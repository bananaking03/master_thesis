module calibration_module #(
    parameter DAC_EXTRA_DATA_WIDTH = 7,
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = 10,
    parameter HISTOGRAM_DATA_WIDTH = 5,
    parameter CALIBRATION_LENGTH = 1800,
    parameter CAL_CONSTANT = 0.1 // scaling constant for histogram difference to threshold update
)(
    input wire clk,
    input wire rst_n,
    input wire [ALGO_DATA_WIDTH-1:0] data_in,
    input wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in,
    output wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_out,
    output wire dither_out
);

    // Internal signals
    localparam EXTRA_WIDTH_FROM_ALGO_TO_DAC = DAC_EXTRA_DATA_WIDTH + DAC_CTRL_WIDTH - ALGO_DATA_WIDTH;
    localparam NUM_HIST = 1 << ALGO_DATA_WIDTH; // number of histogram bins / thresholds
    localparam NUM_THRESHOLDS = NUM_HIST;
    reg [HISTOGRAM_DATA_WIDTH-1:0] hist_plus_data [0:NUM_HIST-1];
    reg [HISTOGRAM_DATA_WIDTH-1:0] hist_minus_data [0:NUM_HIST-1];

    // Generate dither signal
    dither_generator dither_gen (
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
    
    // Accumulate histogram data
    histogram #(
        .DATA_WIDTH(ALGO_DATA_WIDTH)
    ) hist_plus_inst (
        .clk(clk),
        .rst_n(rst_n),
        .acc_en(prev_dither_out),
        .data_in(data_in -1),
        .data_out(hist_plus_data)
    );

    histogram #(
        .DATA_WIDTH(ALGO_DATA_WIDTH)
    ) hist_minus_inst (
        .clk(clk),
        .rst_n(rst_n),
        .acc_en(~prev_dither_out),
        .data_in(data_in),
        .data_out(hist_minus_data)
    );


    // Calculate histogram difference element-wise and multiply with cal_constant
    wire [HISTOGRAM_DATA_WIDTH-1:0] histogram_diff [0:NUM_HIST-1];

    genvar i;
    generate
        for (i = 0; i < NUM_HIST; i = i + 1) begin : gen_hist_diff
            assign histogram_diff[i] = (hist_plus_data[i] - hist_minus_data[i]) * CAL_CONSTANT; // cal_constant is a predefined constant for scaling
        end
    endgenerate

    // Registers to store the thresholds
    reg [HISTOGRAM_DATA_WIDTH-1:0] delta_thresholds_reg [0:NUM_HIST-1];

    // update thresholds after a calibration length of cycles
    // counter width uses a safe 32-bit reg (CALIBRATION_LENGTH fits in 32 bits)
    reg [31:0] calib_cnt;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset delta thresholds to initial values and clear counter
            integer j;
            for (j = 0; j < NUM_THRESHOLDS; j = j + 1) begin
                delta_thresholds_reg[j] <=  2**EXTRA_WIDTH_FROM_ALGO_TO_DAC - 1; // difference between each initial threshold
            end
            calib_cnt <= 32'd0;
        end else begin
            if (calib_cnt == (CALIBRATION_LENGTH - 1)) begin
                // reached calibration length: update thresholds based on histogram_diff
                integer k;
                for (k = 0; k < NUM_HIST; k = k + 1) begin
                    delta_thresholds_reg[k] <= delta_thresholds_reg[k] + histogram_diff[k]; // Update with calculated difference
                end
                calib_cnt <= 32'd0;
            end else begin
                // increment calibration cycle counter
                calib_cnt <= calib_cnt + 1;
            end
        end
    end

    // Select thresholds to generate DAC control output
    // DAC_ctrl_bin = sum(delta_thresholds_reg[0 : floor(DAC_ctrl_in/4)]) + delta_thresholds_reg[ceil(DAC_ctrl_in/4)] / (DAC_ctrl_in % 4)
    
    wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_index_lower = DAC_ctrl_in >> 2;  // floor(DAC_ctrl_in / 4)
    wire [1:0] DAC_ctrl_fraction = DAC_ctrl_in & 3;  // DAC_ctrl_in % 4
    wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_index_upper = DAC_ctrl_index_lower + (|DAC_ctrl_fraction ? 1 : 0);  // ceil(DAC_ctrl_in / 4)
    
    reg [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH - 1:0] DAC_ctrl_bin_sum;
    wire [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH - 1:0] DAC_ctrl_bin_frac;
    
    // Sum from index 0 to floor(DAC_ctrl_in/4)
    always @(*) begin
        integer sum_idx;
        DAC_ctrl_bin_sum = 0;
        for (sum_idx = 0; sum_idx <= DAC_ctrl_index_lower; sum_idx = sum_idx + 1) begin
            DAC_ctrl_bin_sum = DAC_ctrl_bin_sum + delta_thresholds_reg[sum_idx];
        end
    end
    
    // Fractional part: delta_thresholds_reg[ceil_index] / (DAC_ctrl_in % 4)
    assign DAC_ctrl_bin_frac = (DAC_ctrl_fraction == 0) ? 0 : (delta_thresholds_reg[DAC_ctrl_index_upper] / DAC_ctrl_fraction);
    
    wire [HISTOGRAM_DATA_WIDTH + DAC_CTRL_WIDTH - 1:0] DAC_ctrl_bin = DAC_ctrl_bin_sum + DAC_ctrl_bin_frac;

    // Store the final DAC control output in thermometer code format

    // Combinational logic to assign output data
    // Convert DAC_ctrl_bin to thermometer code: bit i = 1 if i < DAC_ctrl_bin, else 0
    genvar therm_idx;
    generate
        for (therm_idx = 0; therm_idx < DAC_CTRL_WIDTH; therm_idx = therm_idx + 1) begin : gen_thermometer
            assign DAC_ctrl_out[therm_idx] = (therm_idx < DAC_ctrl_bin);
        end
    endgenerate

endmodule