module unary_popcount #(
    parameter N = 8191,
    parameter CODE_WIDTH = 13
)(
    input wire [N-1:0] unary_ctrl,
    output reg [CODE_WIDTH-1:0] dac_code
);

    integer i;

    always @(*) begin
        dac_code = 0;

        for (i = 0; i < N; i = i + 1) begin
            dac_code = dac_code + unary_ctrl[i];
        end
    end

endmodule
