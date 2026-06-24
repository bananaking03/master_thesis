`timescale 1ns/1ps

module tb_calibration_module;

    parameter ALGO_DATA_WIDTH = 8;
    parameter DAC_CTRL_WIDTH = 10;
    parameter CLK_PERIOD = 10; // ns
    parameter DAC_EXTRA_DATA_WIDTH = 7; // Added parameter for extra DAC control bits

    reg clk;
    reg rst_n;
    reg [ALGO_DATA_WIDTH-1:0] data_in;
    reg [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in;
    wire [2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)-1:0] DAC_ctrl_out;
    wire dither_out;

    // Instantiate DUT with default parameters
    calibration_module #(
        .DAC_EXTRA_DATA_WIDTH(7),
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .DAC_CTRL_WIDTH(DAC_CTRL_WIDTH),
        .HISTOGRAM_DATA_WIDTH(8),
        .CALIBRATION_LENGTH(1800),
        .CAL_CONSTANT(0.1)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .data_in(data_in),
        .DAC_ctrl_in(DAC_ctrl_in),
        .DAC_ctrl_out(DAC_ctrl_out),
        .dither_out(dither_out)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end

    // Test scenario
    localparam INPUT_LEN = 10000;
    integer sample_idx;
    integer ones_count;
    integer bit_idx;
    real sine_val;

    initial begin
        // Waveform dump
        // $dumpfile("calibration_module_tb.vcd");
        // $dumpvars(0, tb_calibration_module);

        // init
        rst_n = 0;
        data_in = 0;
        DAC_ctrl_in = 0;
        sample_idx = 0;

        // release reset after a few cycles
        #(CLK_PERIOD*5);
        rst_n = 1;

        // Apply sine-wave stimulus to data_in and cycle DAC_ctrl_in through 1..9
        while (sample_idx < INPUT_LEN) begin
            @(posedge clk);
            sine_val = 128.0 + 127.0 * $sin(2.0 * 3.141592653589793 * sample_idx / 100.0);
            data_in = $rtoi(sine_val);
            DAC_ctrl_in = (sample_idx % 9) + 1;

            ones_count = 0;
            for (bit_idx = 0; bit_idx < (2**(DAC_CTRL_WIDTH+DAC_EXTRA_DATA_WIDTH)); bit_idx = bit_idx + 1) begin
                ones_count = ones_count + DAC_ctrl_out[bit_idx];
            end
            $display("%0t: sample=%0d DAC_ctrl_in=%0d DAC_ctrl_out_ones=%0d", $time, sample_idx, DAC_ctrl_in, ones_count);

            sample_idx = sample_idx + 1;
        end

        // Allow a few extra cycles for the DUT to settle
        repeat (10) @(posedge clk);
        $finish;
    end

endmodule
