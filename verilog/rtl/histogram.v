module histogram #(
    parameter ALGO_DATA_WIDTH = 8,
    parameter HISTOGRAM_DATA_WIDTH = 5
) 
(
    input wire clk,
    input wire rst_n,
    input wire acc_en,
    input wire [ALGO_DATA_WIDTH-1:0] data_in,
    output reg [HISTOGRAM_DATA_WIDTH-1:0] data_out [(2**ALGO_DATA_WIDTH)-1:0]
);

    integer i;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset histogram data
            for (i = 0; i < (2**ALGO_DATA_WIDTH); i = i + 1) begin
                data_out[i] <= 0;
            end
        end else if (acc_en) begin
            // Accumulate histogram data based on input data
            data_out[data_in] <= data_out[data_in] + 1; // Increment the count for the corresponding bin
        end
    end


endmodule