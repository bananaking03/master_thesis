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

    reg [ALGO_DATA_WIDTH-1:0] data_mem [0:INPUT_LEN-1];
    reg [DAC_CTRL_WIDTH-1:0] dacctrl_mem [0:INPUT_LEN-1];
    integer sample_idx;

    initial begin
        // Waveform dump
        $dumpfile("calibration_module_tb.vcd");
        $dumpvars(0, tb_calibration_module);

        // Load stimulus from MATLAB-generated memory files
        $readmemh("verilog/tb/calibration_data_in.mem", data_mem);
        $readmemh("verilog/tb/calibration_DAC_ctrl_in.mem", dacctrl_mem);

        // init
        rst_n = 0;
        data_in = 0;
        DAC_ctrl_in = 0;
        sample_idx = 0;

        // release reset after a few cycles
        #(CLK_PERIOD*5);
        rst_n = 1;

        // Apply stimulus from the file
        while (sample_idx < INPUT_LEN) begin
            @(posedge clk);
            data_in = data_mem[sample_idx];
            DAC_ctrl_in = dacctrl_mem[sample_idx];
            sample_idx = sample_idx + 1;
        end

        // Allow a few extra cycles for the DUT to settle
        repeat (10) @(posedge clk);
        $finish;
    end

endmodule
