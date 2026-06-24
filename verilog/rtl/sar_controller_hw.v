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
    output reg  [N-1:0]     cdac_bits,    // drive switches: 1 -> connect corresponding cap to Vref (trial)
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

    integer i;

    // Sequential: state, counters, outputs
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= IDLE;
            cdac_bits <= {N{1'b0}};
            sample <= 1'b0;
            busy <= 1'b0;
            done <= 1'b0;
            bit_idx <= N>0 ? N-1 : 0;
            settle_cnt <= 0;
        end else begin
            state <= next_state;

            case (state)
                IDLE: begin
                    sample <= 1'b0;
                    done <= 1'b0;
                    busy <= 1'b0;
                    if (start) begin
                        busy <= 1'b1;
                        cdac_bits <= {N{1'b0}};
                        bit_idx <= N>0 ? N-1 : 0;
                    end
                end
                SAMPLE: begin
                    sample <= 1'b1;
                end
                CONVERT: begin
                    sample <= 1'b0;
                    // settle counter handled below
                end
                DONE: begin
                    done <= 1'b1;
                    busy <= 1'b0;
                end
            endcase
        end
    end

    // Single-block FSM and conversion control (synchronous)
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= IDLE;
            cdac_bits <= {N{1'b0}};
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
                        cdac_bits <= {N{1'b0}};
                        bit_idx <= N>0 ? N-1 : 0;
                        state <= SAMPLE;
                    end
                end

                SAMPLE: begin
                    // assert sampling for one cycle
                    sample <= 1'b1;
                    // prepare first trial bit and settle counter on next cycle
                    cdac_bits <= {N{1'b0}} | (1 << bit_idx);
                    settle_cnt <= SETTLE_CYCLES;
                    state <= CONVERT;
                end

                CONVERT: begin
                    sample <= 1'b0;
                    if (settle_cnt > 0) begin
                        settle_cnt <= settle_cnt - 1;
                    end else begin
                        // comparator sampled, update the bit
                        if (cmp) begin
                            cdac_bits[bit_idx] <= 1'b1;
                        end else begin
                            cdac_bits[bit_idx] <= 1'b0;
                        end

                        if (bit_idx == 0) begin
                            state <= DONE;
                        end else begin
                            // move to next bit: decrement and start next trial
                            bit_idx <= bit_idx - 1;
                            cdac_bits <= cdac_bits | (1 << (bit_idx - 1));
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

endmodule
