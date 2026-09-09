`timescale 1ns/1ps
// Synthesizable SAR controller that drives CDAC selection bits and sampling phase
module sar_controller_hw #(
    parameter N = 10,
    parameter SETTLE_CYCLES = 2
)(
    input  wire             clk,
    input  wire             rst_n,
    input  wire             start,
    input  wire             cmp,          // comparator result: 1 if Vin >= Vdac
    output reg  [2**N-1:0]     cdac_bits,    // drive switches: 1 -> connect corresponding cap to Vref (trial)
	output wire  [2**N-1:0]  neg_cdac_bits,
    output reg              sample,       // sampling phase (asserted during sample)
    output reg              busy,
    output reg              done
);

    // FSM states
    localparam IDLE  = 2'd0;
    localparam SAMPLE = 2'd1;
    localparam CONVERT = 2'd2;
    localparam DONE  = 2'd3;

    reg [1:0] state, next_state;
    reg [$clog2(N)-1:0] bit_idx;
    reg [$clog2(SETTLE_CYCLES+1)-1:0] settle_cnt;
    reg [N-1:0] cdac_bits_bin;

    integer i, code_val;

    // Single-block FSM and conversion control (synchronous)
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= IDLE;
            cdac_bits_bin <= {N{1'b0}};
            sample <= 1'b0;
            busy <= 1'b0;
            done <= 1'b0;
            bit_idx <= N>0 ? N-1 : 0;
            settle_cnt <= 0;
        end else begin
            case (state)
                IDLE: begin
                    sample <= 1'b0;
                    done <= 1'b0;
                    busy <= 1'b0;
                    if (start) begin
                        busy <= 1'b1;
                        cdac_bits_bin <= {N{1'b0}};
                        bit_idx <= N>0 ? N-1 : 0;
                        state <= SAMPLE;
                    end
                end

                SAMPLE: begin
                    // assert sampling for one cycle
                    sample <= 1'b1;
                    // prepare first trial bit and settle counter on next cycle
                    cdac_bits_bin <= {N{1'b0}} | (1 << bit_idx);
                    settle_cnt <= SETTLE_CYCLES;
                    state <= CONVERT;
                end

                CONVERT: begin
                    sample <= 1'b0;
                    if (settle_cnt > 0) begin
                        settle_cnt <= settle_cnt - 1;
                    end else begin
                        // comparator sampled, update the bit and prepare next trial
                        if (bit_idx == 0) begin
                            // final bit: set or clear LSB according to comparator
                            if (cmp) begin
                                cdac_bits_bin <= cdac_bits_bin | (1 << bit_idx);
                            end else begin
                                cdac_bits_bin <= cdac_bits_bin & ~(1 << bit_idx);
                            end
                            state <= DONE;
                        end else begin
                            // set/clear current bit, and set next trial bit (bit_idx-1)
                            if (cmp) begin
                                cdac_bits_bin <= (cdac_bits_bin | (1 << bit_idx)) | (1 << (bit_idx - 1));
                            end else begin
                                cdac_bits_bin <= (cdac_bits_bin & ~(1 << bit_idx)) | (1 << (bit_idx - 1));
                            end
                            bit_idx <= bit_idx - 1;
                            settle_cnt <= SETTLE_CYCLES;
                            state <= CONVERT;
                        end
                    end
                end

                DONE: begin
                    done <= 1'b1;
                    busy <= 1'b0;
                    // wait for start deassertion to go idle
                    if (!start) begin
                        state <= IDLE;
                        done <= 1'b0;
                    end
                end
            endcase
        end
    end

    // Convert the binary-style SAR decision word into a unary CDAC selection vector.
    always @(*) begin
        cdac_bits = {2**N{1'b0}};
        code_val = 0;

        for (i = 0; i < N; i = i + 1) begin
            if (cdac_bits_bin[i]) begin
                code_val = code_val + (1 << i);
            end
        end

        for (i = 0; i < 2**N; i = i + 1) begin
            if (i < code_val) begin
                cdac_bits[i] = 1'b1;
            end else begin
                cdac_bits[i] = 1'b0;
            end
        end
    end

	assign neg_cdac_bits = ~cdac_bits;

endmodule