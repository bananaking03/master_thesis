`timescale 1ns/1ps

module tb_calibration_module;

    parameter ALGO_DATA_WIDTH = 8;
    parameter DAC_CTRL_WIDTH = 10;
    parameter CLK_PERIOD = 10; // ns

    reg clk;
    reg rst_n;
    reg [ALGO_DATA_WIDTH-1:0] data_in;
    reg [DAC_CTRL_WIDTH-1:0] DAC_ctrl_in;
    wire [DAC_CTRL_WIDTH-1:0] DAC_ctrl_out;
    wire dither_out;

    // Instantiate DUT with default parameters
    calibration_module #(
        .DAC_EXTRA_DATA_WIDTH(7),
        .ALGO_DATA_WIDTH(ALGO_DATA_WIDTH),
        .DAC_CTRL_WIDTH(DAC_CTRL_WIDTH),
        .HISTOGRAM_DATA_WIDTH(5),
        .CALIBRATION_LENGTH(1800)
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
    integer cycle_cnt;
    integer print_cnt;
    initial begin
        // Waveform dump
        $dumpfile("calibration_module_tb.vcd");
        $dumpvars(0, tb_calibration_module);

        // init
        rst_n = 0;
        data_in = 0;
        DAC_ctrl_in = 0;
        cycle_cnt = 0;
        print_cnt = 0;

        // release reset after a few cycles
        #(CLK_PERIOD*5);
        rst_n = 1;

        // Run for a while to cover at least a couple calibration cycles
        while (cycle_cnt < 20000) begin
            @(posedge clk);
            cycle_cnt = cycle_cnt + 1;

            // simple stimulus: ramp data_in and DAC_ctrl_in to exercise logic
            data_in <= data_in + 1;
            DAC_ctrl_in <= DAC_ctrl_in + 1;

            // occasional logging
            print_cnt = print_cnt + 1;
            if (print_cnt == 1000) begin
                $display("time=%0t cycle=%0d data_in=%0d DAC_ctrl_in=%0d dither=%b DAC_ctrl_out=%b", $time, cycle_cnt, data_in, DAC_ctrl_in, dither_out, DAC_ctrl_out);
                print_cnt = 0;
            end
        end

        $display("Testbench finished at time %0t", $time);
        $finish;
    end

endmodule
