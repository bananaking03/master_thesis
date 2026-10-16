module digi_top #(
    parameter N = 10,
    parameter ALGO_DATA_WIDTH = 8,
    parameter DAC_CTRL_WIDTH = 8,
    parameter HISTOGRAM_DATA_WIDTH = 5
)(
    input wire clk,
    input wire rst_n,
    input wire start,
    input wire cmp, // comparator result: 1 if Vin >= Vdac
    output wire dither_out,
    output reg [2**(N-1)-1:0] cdac_bits_ctrl, // drive switches: 1 -> connect corresponding cap to Vref (trial)
    output reg [2**(N-1)-1:0] neg_cdac_bits_ctrl,
    output reg sample, // sampling phase (asserted during sample)
);

    reg [2**(N-1)-1:0] cdac_bits_bin;

    // Instantiate SAR logic
    SAR_logic #(
        .N(N)
    ) sar_logic_inst (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .cmp(cmp),
        .cdac_bits(cdac_bits),
        .neg_cdac_bits(neg_cdac_bits),
        .sample(sample),
        .busy(busy),
        .done(done),
        .cdac_bits_bin(cdac_bits_bin)
    );

    // Instantiate calibration module
    calibration_module #(
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .DAC_CTRL_WIDTH(DAC_CTRL_WIDTH),
        .HISTOGRAM_DATA_WIDTH(HISTOGRAM_DATA_WIDTH)
    ) calibration_module_inst (
        .clk(clk),
        .rst_n(rst_n),
        .data_in(cdac_bits_bin),
        .DAC_ctrl_out(cdac_bits_ctrl),
        .DAC_ctrl_out_neg(neg_cdac_bits_ctrl),
        .dither_out(dither_out)
    );
endmodule