module unary_popcount #(
    parameter int N = 8192,
    parameter int CODE_WIDTH = 13
) (
    input  logic [N-1:0] unary_ctrl,
    output logic [CODE_WIDTH-1:0] dac_code
);

    integer i;
    integer count;

    always_comb begin
        count = 0;

        for (i = 0; i < N; i = i + 1)
            count = count + unary_ctrl[i];

        dac_code = count;
    end

endmodule