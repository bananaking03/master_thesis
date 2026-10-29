module calibration_module #(
    parameter DAC_EXTRA_DATA_WIDTH = 7,
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = 10,
    parameter HISTOGRAM_DATA_WIDTH = 8,
    parameter CALIBRATION_LENGTH = 1800,
    parameter CAL_CONSTANT = 0.1, // scaling constant for histogram difference to threshold update (kept for compatibility)
    parameter integer CAL_CONSTANT_NUM = 1, // integer numerator for fixed-point scaling (CAL_CONSTANT_NUM/CAL_CONSTANT_DEN)
    parameter integer CAL_CONSTANT_DEN = 10, // integer denominator for fixed-point scaling (default 1/10 = 0.1)
    parameter integer LAMBDA_REG_NUM = 1, // lambda_reg numerator (default 1/10 = 0.1)
    parameter integer LAMBDA_REG_DEN = 10,
    parameter integer SOLVER_FRAC_BITS = 16 // fractional bits used by the matrix solver
)(
    input wire clk,
    input wire sar_done,
    input wire sample,
    input wire rst_n,
    input wire signed [ALGO_DATA_WIDTH-1:0] data_in,
    input wire signed [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in,
    output wire [(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out_bin,
    output wire [(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out_neg_bin,
    output wire dither_out
);

    // Internal signals
    localparam EXTRA_WIDTH_FROM_ALGO_TO_DAC = DAC_EXTRA_DATA_WIDTH + DAC_CTRL_WIDTH - ALGO_DATA_WIDTH;
    localparam NUM_HIST = 1 << ALGO_DATA_WIDTH; // number of histogram bins / thresholds
    localparam NUM_THRESHOLDS = NUM_HIST;
    localparam HISTOGRAM_ROTATION = (NUM_HIST / 2 + 1) % NUM_HIST;
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_plus_data [0:NUM_HIST-1];
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_minus_data [0:NUM_HIST-1];
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_plus_raw [0:NUM_HIST-1];
    wire [HISTOGRAM_DATA_WIDTH-1:0] hist_minus_raw [0:NUM_HIST-1];

    genvar hist_map_idx;
    generate
        for (hist_map_idx = 0; hist_map_idx < NUM_HIST; hist_map_idx = hist_map_idx + 1) begin : gen_histogram_reorder
            assign hist_plus_data[hist_map_idx] = hist_plus_raw[(hist_map_idx + HISTOGRAM_ROTATION) % NUM_HIST];
            assign hist_minus_data[hist_map_idx] = hist_minus_raw[(hist_map_idx + HISTOGRAM_ROTATION) % NUM_HIST];
        end
    endgenerate

    // Generate dither signal
    dither_gen dither_gen (
        .clk(sar_done),
        .rst_n(rst_n),
        .dither_out(dither_out)
    );
    
    // Register previous-cycle dither value for use by other blocks
    reg prev_dither_out;
    always @(posedge sar_done or negedge rst_n) begin
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
        .clk(sar_done),
        .rst_n(rst_n),
        .acc_en(dither_out),
        .reset_bins(thresholds_update_now),
        .data_in(data_in_minus1),
        .data_out(hist_plus_raw)
    );

    histogram #(
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .HISTOGRAM_DATA_WIDTH(HISTOGRAM_DATA_WIDTH)
    ) hist_minus_inst (
        .clk(sar_done),
        .rst_n(rst_n),
        .acc_en(~dither_out),
        .reset_bins(thresholds_update_now),
        .data_in(data_in),
        .data_out(hist_minus_raw)
    );


    // Calculate histogram difference element-wise (signed)
    wire signed [HISTOGRAM_DATA_WIDTH:0] histogram_diff [0:NUM_HIST-1];

    genvar i;
    generate
        for (i = 0; i < NUM_HIST; i = i + 1) begin : gen_hist_diff
            // Keep the raw difference for matching; scale the matched result afterward.
            assign histogram_diff[i] = $signed(hist_plus_data[i]) - $signed(hist_minus_data[i]);
        end
    endgenerate

    // MATLAB solves (error_matrix + lambda_reg*I) \ H_delta.
    // The matrix is tridiagonal, with 0.5 + lambda_reg on the diagonal and
    // -0.25 on both neighboring diagonals. Multiplying the matrix by
    // 4*LAMBDA_REG_DEN gives integer coefficients without changing the solution.
    localparam integer MATRIX_DIAG = (2 * LAMBDA_REG_DEN) + (4 * LAMBDA_REG_NUM);
    localparam integer MATRIX_OFF_DIAG = -LAMBDA_REG_DEN;

    // Fixed-point arrays used by the Thomas algorithm for A*x = H_delta.
    reg signed [63:0] solver_cprime [0:NUM_HIST-1];
    reg signed [63:0] solver_dprime [0:NUM_HIST-1];
    reg signed [63:0] matched_histogram [0:NUM_HIST-1];
    wire signed [63:0] scaled_matched_histogram [0:NUM_HIST-1];
    wire signed [63:0] matched_histogram_update [0:NUM_HIST-1];
    integer solver_index;
    reg signed [63:0] solver_denom;
    reg signed [63:0] solver_rhs;

    genvar scale_idx;
    generate
        for (scale_idx = 0; scale_idx < NUM_HIST; scale_idx = scale_idx + 1) begin : gen_scaled_matched_histogram
            assign scaled_matched_histogram[scale_idx] =
                (((matched_histogram[scale_idx] * CAL_CONSTANT_NUM) / CAL_CONSTANT_DEN)/CALIBRATION_LENGTH) <<< (ALGO_DATA_WIDTH - 1);
            assign matched_histogram_update[scale_idx] =
                scaled_matched_histogram[scale_idx] >>> SOLVER_FRAC_BITS;
        end
    endgenerate

    always @(*) begin
        // Forward elimination for the tridiagonal error matrix.
        solver_cprime[0] = (MATRIX_OFF_DIAG <<< SOLVER_FRAC_BITS) / MATRIX_DIAG;
        solver_rhs = $signed(histogram_diff[0]) <<< SOLVER_FRAC_BITS;
        solver_dprime[0] = (solver_rhs <<< SOLVER_FRAC_BITS) / MATRIX_DIAG;

        for (solver_index = 1; solver_index < NUM_HIST; solver_index = solver_index + 1) begin
            solver_denom = MATRIX_DIAG -
                           (MATRIX_OFF_DIAG * solver_cprime[solver_index - 1] >>> SOLVER_FRAC_BITS);
            solver_cprime[solver_index] =
                (MATRIX_OFF_DIAG <<< SOLVER_FRAC_BITS) / solver_denom;
            solver_rhs = ($signed(histogram_diff[solver_index]) <<< SOLVER_FRAC_BITS) -
                         (MATRIX_OFF_DIAG * solver_dprime[solver_index - 1] >>> SOLVER_FRAC_BITS);
            solver_dprime[solver_index] =
                (solver_rhs <<< SOLVER_FRAC_BITS) / solver_denom;
        end

        // Back substitution produces the regularized histogram difference.
        matched_histogram[NUM_HIST - 1] = solver_dprime[NUM_HIST - 1];
        for (solver_index = NUM_HIST - 2; solver_index >= 0; solver_index = solver_index - 1) begin
            matched_histogram[solver_index] = solver_dprime[solver_index] -
                (solver_cprime[solver_index] * matched_histogram[solver_index + 1] >>> SOLVER_FRAC_BITS);
        end
    end

    // Registers to store the thresholds
    reg [EXTRA_WIDTH_FROM_ALGO_TO_DAC+1:0] delta_thresholds_reg [0:NUM_HIST-1];

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset delta thresholds to initial values and clear counter
            integer j;
            for (j = 0; j < NUM_THRESHOLDS; j = j + 1) begin
                delta_thresholds_reg[j] <= 2**EXTRA_WIDTH_FROM_ALGO_TO_DAC;
            end
            calib_cnt <= 32'd0;
        end else begin
            if (thresholds_update_now) begin
                // reached calibration length: update thresholds based on histogram_diff
                integer k;
                signed [EXTRA_WIDTH_FROM_ALGO_TO_DAC+1:0] signed_result;
                for (k = 0; k < NUM_HIST; k = k + 1) begin
                    // Convert the fixed-point solver result back to an integer.
                    signed_result = $signed({1'b0, delta_thresholds_reg[k]}) +
                                    matched_histogram_update[k];
                    if (signed_result < 0) // prevent overflow and underflow of thresholds
                        delta_thresholds_reg[k] <= {EXTRA_WIDTH_FROM_ALGO_TO_DAC+2{1'b0}};
                    else
                        delta_thresholds_reg[k] <= signed_result;
                end
                calib_cnt <= 32'd0;
            end else if (sar_done) begin
                // increment calibration cycle counter
                calib_cnt <= calib_cnt + 1;
            end
        end
    end

    // Select thresholds to generate DAC control output
    // DAC_ctrl_in is signed; negative codes use the negative rail and retain
    // the controller's one-step offset when converting to a magnitude.
    // The magnitude is partitioned in groups of 4 steps:
    // 0 -> 0
    // 1 -> delta[0]/4, 2 -> delta[0]*2/4, 3 -> delta[0]*3/4, 4 -> delta[0]
    // 5 -> delta[0] + delta[1]/4, 6 -> delta[0] + delta[1]*2/4, etc.
    wire signed [DAC_CTRL_WIDTH:0] DAC_ctrl_extended = {DAC_ctrl_in[DAC_CTRL_WIDTH-1], DAC_ctrl_in};
    wire [DAC_CTRL_WIDTH:0] DAC_ctrl_magnitude = DAC_ctrl_in[DAC_CTRL_WIDTH-1]
        ? -DAC_ctrl_extended
        : DAC_ctrl_extended;
    wire [DAC_CTRL_WIDTH:0] DAC_ctrl_selection_code = DAC_ctrl_magnitude + DAC_ctrl_in[DAC_CTRL_WIDTH-1];
    wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_index_lower = DAC_ctrl_selection_code >> 2;
    wire [1:0] DAC_ctrl_fraction = DAC_ctrl_selection_code[1:0];
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

    // Convert the calibrated threshold count to active-low thermometer code.
    // Negative codes use the negative rail; positive codes use the primary rail.
    wire [2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out,
    wire [2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out_neg,
    genvar therm_idx;
    generate
        for (therm_idx = 0; therm_idx < (2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)); therm_idx = therm_idx + 1) begin : gen_thermometer
            assign DAC_ctrl_out[therm_idx] = sample || !DAC_ctrl_in[DAC_CTRL_WIDTH-1] || (therm_idx >= DAC_ctrl_bin);
            assign DAC_ctrl_out_neg[therm_idx] = sample || DAC_ctrl_in[DAC_CTRL_WIDTH-1] || (therm_idx >= DAC_ctrl_bin);
        end
    endgenerate

    unary_popcount #(
        .N(2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)),
        .CODE_WIDTH(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)
    ) unary_popcount_inst (
        .in(DAC_ctrl_out),
        .out(DAC_ctrl_out_bin)
    );

    unary_popcount #(
        .N(2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)),
        .CODE_WIDTH(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)
    ) unary_popcount_neg_inst (
        .in(DAC_ctrl_out_neg),
        .out(DAC_ctrl_out_neg_bin)
    );

endmodule