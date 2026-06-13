module histogram #(
    parameter ALGO_DATA_WIDTH = 8,
    parameter HISTOGRAM_DATA_WIDTH = 5,
    parameter NUM_BINS = (1 << ALGO_DATA_WIDTH)
) 
(
    input  logic clk,
    input  logic rst_n,
    input  logic acc_en,
    input  logic [ALGO_DATA_WIDTH-1:0] data_in,
    output logic [HISTOGRAM_DATA_WIDTH-1:0] data_out [0:NUM_BINS-1]
);

    // Internal storage driven by procedural logic
    logic [HISTOGRAM_DATA_WIDTH-1:0] bin_reg [0:NUM_BINS-1];
    integer i;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset histogram data
            for (i = 0; i < NUM_BINS; i = i + 1) begin
                bin_reg[i] <= '{default: '0};
            end
        end else if (acc_en) begin
            // Accumulate histogram data based on input data
            bin_reg[data_in] <= bin_reg[data_in] + 1; // Increment the count for the corresponding bin
        end
    end

    // Drive output ports from internal regs (packed-to-unpacked assignment)
    genvar gi;
    generate
        for (gi = 0; gi < NUM_BINS; gi = gi + 1) begin : gen_drive
            assign data_out[gi] = bin_reg[gi];
        end
    endgenerate

endmodule
