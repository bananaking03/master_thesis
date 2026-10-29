`timescale 1ns/1ps
// Synthesizable SAR controller that drives CDAC selection bits and sampling phase
module SAR_logic #(
    parameter N = 10
)(
    input  wire             clk,
    input  wire             rst_n,
    input  wire             start,
    input  wire             cmp,          // comparator result: 1 if Vin >= Vdac
    output reg              sample,       // sampling phase (asserted during sample)
    output reg              busy,
    output reg              done,
    output reg signed [N-1:0] cdac_bits_bin
);

    // FSM states
    localparam IDLE  = 2'd0;
    localparam SAMPLE = 2'd1;
    localparam CONVERT = 2'd2;
    localparam DONE  = 2'd3;

    reg [1:0] state, next_state;
    reg [$clog2(N)-1:0] bit_idx;
	reg signed [N-1:0] result;

    // Single-block FSM and conversion control (synchronous)
    always @(negedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= IDLE;
            cdac_bits_bin <= {N{1'b0}};
            sample <= 1'b0;
            busy <= 1'b0;
            done <= 1'b0;
            bit_idx <= N>0 ? N-1 : 0;
			result <= {N{1'b0}};
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
                        sample <= 1'b1;
                        state <= SAMPLE;
                    end
                end

                SAMPLE: begin
                    // assert sampling for one cycle
                    sample <= 1'b0;
                    //cdac_bits_bin <= {N{1'b0}} | (1 << bit_idx);
                    cdac_bits_bin <= {N{1'b0}};
					done <= 1'b0;
                    result <= cdac_bits_bin;
                    state <= CONVERT;
                end

                CONVERT: begin
                    // comparator sampled, update the bit and prepare next trial
                    if (bit_idx == 0) begin
                        // final bit: set or clear LSB according to comparator
                        if (cmp) begin
                            cdac_bits_bin <= cdac_bits_bin | (1 << bit_idx);
                        end else begin
                            cdac_bits_bin <= cdac_bits_bin & ~(1 << bit_idx);
                        end
                        done <= 1'b1;
                        if (start) begin
                            state <= SAMPLE;
                            bit_idx <= N>0 ? N-1 : 0;
                        sample <= 1'b1;
                            busy <= 1'b1;
                        end else begin
                            state <= DONE;
                        end
                    end else if (bit_idx == N-1) begin
                        // first bit: set or clear MSB according to comparator, first bit indicates the sign and set next trial bit (bit_idx-1)
                        if (cmp) begin
                            cdac_bits_bin <= cdac_bits_bin & ~(1 << bit_idx) | (1 << (bit_idx - 1));
                        end else begin
                            cdac_bits_bin <= cdac_bits_bin | (1 << bit_idx) | (1 << (bit_idx - 1));
                        end
                        bit_idx <= bit_idx - 1;
                        state <= CONVERT;
                    end else begin
                        // set/clear current bit, and set next trial bit (bit_idx-1)
                        if (cmp) begin
                            cdac_bits_bin <= (cdac_bits_bin | (1 << bit_idx)) | (1 << (bit_idx - 1));
                        end else begin
                            cdac_bits_bin <= (cdac_bits_bin & ~(1 << bit_idx)) | (1 << (bit_idx - 1));
                        end
                        bit_idx <= bit_idx - 1;
                        state <= CONVERT;
                    end
                end

                DONE: begin
                    done <= 1'b1;
                    busy <= 1'b0;
                    // wait for start deassertion to go idle
                    if (!start) begin
                        state <= IDLE;
                        done <= 1'b0;
                    end else begin
					state <= SAMPLE;
                    sample <= 1'b1;
					done <= 1'b0;
					cdac_bits_bin <= {N{1'b0}};
                    bit_idx <= N>0 ? N-1 : 0;
					busy <= 1'b1;
					end
                end
            endcase
        end
    end

endmodule
