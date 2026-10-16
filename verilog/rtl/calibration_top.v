module calibration_module #(
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = 8,
    parameter HISTOGRAM_DATA_WIDTH = 5
)(
    input wire clk,
    input wire rst_n,
    input wire [ALGO_DATA_WIDTH-1:0] data_in,
    output wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_out,
    output wire dither_out
);

    // Internal signals
    reg [(2**ALGO_DATA_WIDTH)-1:0] [HISTOGRAM_DATA_WIDTH-1:0] hist_plus_data;
    reg [(2**ALGO_DATA_WIDTH)-1:0] [HISTOGRAM_DATA_WIDTH-1:0] hist_minus_data;

    // Generate dither signal
    dither_generator dither_gen (
        .clk(clk),
        .rst_n(rst_n),
        .dither_out(dither_out)
    );
    
    // Register previous-cycle dither value for other blocks
    reg prev_dither_out;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            prev_dither_out <= 1'b0;
        else
            prev_dither_out <= dither_out;
    end
    
    // Accumulate histogram data
    histogram_plus #(
        .DATA_WIDTH(ALGO_DATA_WIDTH)
    ) histogram_inst (
        .clk(clk),
        .rst_n(rst_n),
        .dither_in(prev_dither_out),
        .data_in(data_in -1),
        .data_out(hist_plus_data)
    );

    histogram_minus #(
        .DATA_WIDTH(ALGO_DATA_WIDTH)
    ) histogram_inst (
        .clk(clk),
        .rst_n(rst_n),
        .dither_in(prev_dither_out),
        .data_in(data_in),
        .data_out(hist_minus_data)
    );

    // Calculate histogram difference element-wise and multiply with cal_constant
    wire [(2**ALGO_DATA_WIDTH)-1:0] [HISTOGRAM_DATA_WIDTH-1:0] histogram_diff;

    genvar i;
    generate
        for (i = 0; i < (2**ALGO_DATA_WIDTH); i = i + 1) begin : gen_hist_diff
            assign histogram_diff[i] = (hist_plus_data[i] - hist_minus_data[i]) * cal_constant; // cal_constant is a predefined constant for scaling
        end
    endgenerate

    // Combinational logic to assign output data
    assign DAC_ctrl_out = data_reg; // Output the stored data

endmodule