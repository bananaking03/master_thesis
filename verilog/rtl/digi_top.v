module digi_top #(
    parameter N = 10,
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = N,
    parameter DAC_EXTRA_DATA_WIDTH = 4,
    parameter HISTOGRAM_DATA_WIDTH = 5
)(
    input wire clk,
    input wire rst_n,
    input wire start,
    input wire cmp,
    output wire dither_out,
    output wire [16383:0] cdac_bits_ctrl,
    output wire [16383:0] neg_cdac_bits_ctrl,
    output wire sample
);

    wire busy;
    wire sar_done;
    wire signed [N-1:0] cdac_bits_bin;
    wire signed [ALGO_DATA_WIDTH-1:0] calibration_data = cdac_bits_bin;
    wire signed [DAC_CTRL_WIDTH-1:0] calibration_code = cdac_bits_bin;

    SAR_logic #(
        .N(N)
    ) sar_logic_inst (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .cmp(cmp),
        .sample(sample),
        .busy(busy),
        .done(sar_done),
        .cdac_bits_bin(cdac_bits_bin)
    );

    calibration_module #(
        .DAC_EXTRA_DATA_WIDTH(DAC_EXTRA_DATA_WIDTH),
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .DAC_CTRL_WIDTH(DAC_CTRL_WIDTH),
        .HISTOGRAM_DATA_WIDTH(HISTOGRAM_DATA_WIDTH)
    ) calibration_module_inst (
        .clk(clk),
        .sar_done(sar_done),
        .sample(sample),
        .rst_n(rst_n),
        .data_in(calibration_data),
        .DAC_ctrl_in(calibration_code),
        .DAC_ctrl_out(cdac_bits_ctrl),
        .DAC_ctrl_out_neg(neg_cdac_bits_ctrl),
        .dither_out(dither_out)
    );
endmodule