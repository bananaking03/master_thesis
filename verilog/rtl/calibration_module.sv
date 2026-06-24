module calibration_module #(
    parameter DAC_EXTRA_DATA_WIDTH = 7,
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = 10,
    parameter HISTOGRAM_DATA_WIDTH = 8,
    parameter CALIBRATION_LENGTH = 1800,
    parameter CAL_CONSTANT = 0.1, // scaling constant for histogram difference to threshold update (kept for compatibility)
    parameter integer CAL_CONSTANT_NUM = 1, // integer numerator for fixed-point scaling (CAL_CONSTANT_NUM/CAL_CONSTANT_DEN)
    parameter integer CAL_CONSTANT_DEN = 10 // integer denominator for fixed-point scaling (default 1/10 = 0.1)
)(
    input wire clk,
    input wire rst_n,
    input wire [ALGO_DATA_WIDTH-1:0] data_in,
    input wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in,
    output wire [2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out,
    output wire dither_out
);

    // Internal signals
    localparam EXTRA_WIDTH_FROM_ALGO_TO_DAC = DAC_EXTRA_DATA_WIDTH + DAC_CTRL_WIDTH - ALGO_DATA_WIDTH;
    localparam NUM_HIST = 1 << ALGO_DATA_WIDTH; // number of histogram bins / thresholds
    localparam NUM_THRESHOLDS = NUM_HIST;
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_plus_data [0:NUM_HIST-1];
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_minus_data [0:NUM_HIST-1];

    // Generate dither signal
    dither_gen dither_gen (
        .clk(clk),
        .rst_n(rst_n),
        .dither_out(dither_out)
    );
    
    // Register previous-cycle dither value for use by other blocks
    reg prev_dither_out;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            prev_dither_out <= 1'b0;
        else
            prev_dither_out <= dither_out;
    end

    // update thresholds after a calibration length of cycles
    // counter width uses a safe 32-bit reg (CALIBRATION_LENGTH fits in 32 bits)
    reg [31:0] calib_cnt;

    // combinational signal asserted in the same cycle when calib_cnt == CALIBRATION_LENGTH-1
    wire thresholds_update_now = (calib_cnt == (CALIBRATION_LENGTH - 1));

    
    // Accumulate histogram data
    wire [ALGO_DATA_WIDTH-1:0] data_in_minus1;
    assign data_in_minus1 = data_in - 1;

    histogram #(
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .HISTOGRAM_DATA_WIDTH(HISTOGRAM_DATA_WIDTH)
    ) hist_plus_inst (
        .clk(clk),
        .rst_n(rst_n),
        .acc_en(prev_dither_out),
        .reset_bins(thresholds_update_now),
        .data_in(data_in_minus1),
        .data_out(hist_plus_data)
    );

    histogram #(
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .HISTOGRAM_DATA_WIDTH(HISTOGRAM_DATA_WIDTH)
    ) hist_minus_inst (
        .clk(clk),
        .rst_n(rst_n),
        .acc_en(~prev_dither_out),
        .reset_bins(thresholds_update_now),
        .data_in(data_in),
        .data_out(hist_minus_data)
    );


    // Calculate histogram difference element-wise (signed)
    wire signed [HISTOGRAM_DATA_WIDTH:0] histogram_diff [0:NUM_HIST-1];

    genvar i;
    generate
        for (i = 0; i < NUM_HIST; i = i + 1) begin : gen_hist_diff
            // compute signed difference and apply integer fixed-point scaling (NUM/DEN)
            assign histogram_diff[i] = ($signed($signed(hist_plus_data[i]) - $signed(hist_minus_data[i])) * CAL_CONSTANT_NUM) / CAL_CONSTANT_DEN;
        end
    endgenerate

    // Registers to store the thresholds
    reg [EXTRA_WIDTH_FROM_ALGO_TO_DAC+1:0] delta_thresholds_reg [0:NUM_HIST-1];

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset delta thresholds to initial values and clear counter
            integer j;
            for (j = 0; j < NUM_THRESHOLDS; j = j + 1) begin
                delta_thresholds_reg[j] <= 2**EXTRA_WIDTH_FROM_ALGO_TO_DAC - 1;
            end
            calib_cnt <= 32'd0;
        end else begin
            if (thresholds_update_now) begin
                // reached calibration length: update thresholds based on histogram_diff
                integer k;
                signed [EXTRA_WIDTH_FROM_ALGO_TO_DAC+1:0] signed_result;
                for (k = 0; k < NUM_HIST; k = k + 1) begin
                    signed_result = $signed({1'b0, delta_thresholds_reg[k]}) + histogram_diff[k];
                    if (signed_result < 0) // prevent overflow and underflow of thresholds
                        delta_thresholds_reg[k] <= {EXTRA_WIDTH_FROM_ALGO_TO_DAC+2{1'b0}};
                    else
                        delta_thresholds_reg[k] <= signed_result;
                end
                calib_cnt <= 32'd0;
            end else begin
                // increment calibration cycle counter
                calib_cnt <= calib_cnt + 1;
            end
        end
    end

    // Select thresholds to generate DAC control output
    // DAC_ctrl_in is partitioned in groups of 4 steps:
    // 0 -> 0
    // 1 -> delta[0]/4, 2 -> delta[0]*2/4, 3 -> delta[0]*3/4, 4 -> delta[0]
    // 5 -> delta[0] + delta[1]/4, 6 -> delta[0] + delta[1]*2/4, etc.
    wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_index_lower = DAC_ctrl_in >> 2;  // floor(DAC_ctrl_in / 4)
    wire [1:0] DAC_ctrl_fraction = DAC_ctrl_in & 3;  // DAC_ctrl_in % 4
    reg [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH - 1:0] DAC_ctrl_bin;

    always @(*) begin
        reg [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH - 1:0] base_value;
        reg signed [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH:0] slope;
        reg signed [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH:0] frac_value;
        integer i;

        // base_value is sum of all previous delta thresholds
        base_value = 0;
        for (i = 0; i < DAC_ctrl_index_lower && i < NUM_HIST; i = i + 1) begin
            base_value = base_value + delta_thresholds_reg[i];
        end

        if (DAC_ctrl_index_lower < NUM_HIST) begin
            slope = $signed({1'b0, delta_thresholds_reg[DAC_ctrl_index_lower]});
        end else begin
            slope = 0;
        end

        if (DAC_ctrl_fraction == 0) begin
            DAC_ctrl_bin = base_value;
        end else begin
            frac_value = slope * DAC_ctrl_fraction;
            DAC_ctrl_bin = base_value + (frac_value >>> 2);
        end
    end

    // Store the final DAC control output in thermometer code format

    // Combinational logic to assign output data
    // Convert DAC_ctrl_bin to thermometer code: bit i = 1 if i < DAC_ctrl_bin, else 0
    genvar therm_idx;
    generate
        for (therm_idx = 0; therm_idx < (2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)); therm_idx = therm_idx + 1) begin : gen_thermometer
            assign DAC_ctrl_out[therm_idx] = (therm_idx < DAC_ctrl_bin);
        end
    endgenerate

endmodule