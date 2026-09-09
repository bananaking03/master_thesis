`timescale 1ns/1ps
module sar_controller_hw_tb;
    parameter N = 10;

    reg clk = 0;
    reg rst_n = 0;
    reg start = 0;
    reg [N-1:0] vin = 0;

    wire [N-1:0] cdac_bits;
    wire sample;
    wire busy;
    wire done;
    wire cmp;

    // Instantiate controller from rtl folder
    sar_controller_hw #(.N(N), .SETTLE_CYCLES(2)) dut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .cmp(cmp),
        .cdac_bits(cdac_bits),
        .sample(sample),
        .busy(busy),
        .done(done)
    );

    // Behavioral comparator: decide based on vin vs cdac_bits
    assign cmp = (vin >= cdac_bits) ? 1'b1 : 1'b0;

    always #5 clk = ~clk;

    integer i;

    initial begin
        $display("Starting controller TB");
        // $dumpfile("sar_controller_hw_tb.vcd");
        // $dumpvars(0, sar_controller_hw_tb);

        // reset
        rst_n = 0; start = 0; vin = 0;
        #20 rst_n = 1;

        // Apply several random inputs
        for (i = 0; i < 12; i = i + 1) begin
            vin = $urandom_range(0, (1<<N)-1);
            @(posedge clk);
            start = 1;
            @(posedge clk);
            start = 0;
            wait(done == 1);
            $display("vin=%0d cdac_bits=%b", vin, cdac_bits);
            @(posedge clk);
            #2;
        end

        $display("Finished");
        #20 $finish;
    end

endmodule
