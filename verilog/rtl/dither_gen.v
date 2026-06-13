module digher_gen (
    input wire clk,
    input wire rst_n,
    output reg dither_out
);

    reg [15:0] lfsr; // 16-bit LFSR for pseudo-random generation

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            lfsr <= 16'hACE1; // Non-zero seed value
            dither_out <= 0;
        end else begin
            // LFSR feedback taps (example: taps at bits 16, 14, 13, and 11)
            lfsr <= {lfsr[14:0], lfsr[15] ^ lfsr[13] ^ lfsr[12] ^ lfsr[10]};
            dither_out <= lfsr[0]; // Output the least significant bit as dither signal
        end
    end
endmodule