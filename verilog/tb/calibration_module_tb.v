`timescale 1ns/1ps

module tb_calibration_module;

    localparam ALGO_DATA_WIDTH = 4;
    localparam DAC_CTRL_WIDTH = 4;
    localparam DAC_EXTRA_DATA_WIDTH = 2;
    localparam OUTPUT_WIDTH = 2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH);

    reg clk = 1'b0;
    reg rst_n = 1'b0;
    reg sar_done = 1'b0;
    reg sample = 1'b0;
    reg signed [ALGO_DATA_WIDTH-1:0] data_in = 0;
    reg signed [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in = 0;
    wire [OUTPUT_WIDTH-1:0] DAC_ctrl_out;
    wire [OUTPUT_WIDTH-1:0] DAC_ctrl_out_neg;
    wire dither_out;

    calibration_module #(
        .DAC_EXTRA_DATA_WIDTH(DAC_EXTRA_DATA_WIDTH),
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .DAC_CTRL_WIDTH(DAC_CTRL_WIDTH),
        .HISTOGRAM_DATA_WIDTH(4)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .sar_done(sar_done),
        .sample(sample),
        .data_in(data_in),
        .DAC_ctrl_in(DAC_ctrl_in),
        .DAC_ctrl_out(DAC_ctrl_out),
        .DAC_ctrl_out_neg(DAC_ctrl_out_neg),
        .dither_out(dither_out)
    );

    always #5 clk = ~clk;

    function integer count_zero_bits;
        input [OUTPUT_WIDTH-1:0] value;
        integer bit_index;
        begin
            count_zero_bits = 0;
            for (bit_index = 0; bit_index < OUTPUT_WIDTH; bit_index = bit_index + 1)
                if (!value[bit_index])
                    count_zero_bits = count_zero_bits + 1;
        end
    endfunction

    initial begin
        #12 rst_n = 1'b1;
        #1;
        if (count_zero_bits(DAC_ctrl_out) != 0 || count_zero_bits(DAC_ctrl_out_neg) != 0)
            $fatal(1, "zero input must leave both rails inactive");

        sample = 1'b1;
        DAC_ctrl_in = -4'sd1;
        #1;
        if (count_zero_bits(DAC_ctrl_out) != 0 || count_zero_bits(DAC_ctrl_out_neg) != 0)
            $fatal(1, "sample phase must leave both rails inactive");

        sample = 1'b0;
        DAC_ctrl_in = 4'sd4;
        #1;
        if (count_zero_bits(DAC_ctrl_out) != 3 || count_zero_bits(DAC_ctrl_out_neg) != 0)
            $fatal(1, "positive code must select the primary output rail");

        DAC_ctrl_in = -4'sd1;
        #1;
        if (count_zero_bits(DAC_ctrl_out) != 0 || count_zero_bits(DAC_ctrl_out_neg) != 1)
            $fatal(1, "negative code must select the negative output rail with its offset");

        DAC_ctrl_in = -4'sd8;
        #1;
        if (count_zero_bits(DAC_ctrl_out) != 0 || count_zero_bits(DAC_ctrl_out_neg) != 6)
            $fatal(1, "minimum signed code magnitude was not mapped correctly");

        $display("Signed calibration rail checks passed");
        $finish;
    end

endmodule
