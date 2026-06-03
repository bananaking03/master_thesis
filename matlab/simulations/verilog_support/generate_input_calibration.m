function [data_in, dither_bits, DAC_ctrl_in, thresholds, analog_in] = generate_input_calibration(cal_cycles, N_bits, DAC_CTRL_WIDTH, Vhigh, Vlow, mismatch_LSC, LSC, fs, saveFile)
%GENERATE_INPUT_CALIBRATION Create ADC input and dither for calibration.
%   [data_in, dither_bits, DAC_ctrl_in, thresholds, analog_in] =
%       generate_input_calibration(cal_cycles, N_bits, DAC_CTRL_WIDTH,
%       Vhigh, Vlow, mismatch_LSC, LSC, fs, saveFile)
%
%   Uses existing helper functions from matlab/help_funcs to generate
%   capacitor mismatch and threshold values, then converts a sine-wave
%   input into digital ADC codes with dither.
%
%   Optional arguments are:
%       cal_cycles     default 10000
%       N_bits         default 8
%       DAC_CTRL_WIDTH default N_bits + 2
%       Vhigh          default 1
%       Vlow           default -1
%       mismatch_LSC   default 0.02
%       LSC            default 1
%       fs             default 48000
%       saveFile       default true

    if nargin < 1 || isempty(cal_cycles)
        cal_cycles = 10000;
    end
    if nargin < 2 || isempty(N_bits)
        N_bits = 8;
    end
    if nargin < 3 || isempty(DAC_CTRL_WIDTH)
        DAC_CTRL_WIDTH = N_bits + 2;
    end
    if nargin < 4 || isempty(Vhigh)
        Vhigh = 1;
    end
    if nargin < 5 || isempty(Vlow)
        Vlow = -1;
    end
    if nargin < 6 || isempty(mismatch_LSC)
        mismatch_LSC = 0.02;
    end
    if nargin < 7 || isempty(LSC)
        LSC = 1;
    end
    if nargin < 8 || isempty(fs)
        fs = 48000;
    end
    if nargin < 9 || isempty(saveFile)
        saveFile = true;
    end

    scriptDir = fileparts(mfilename('fullpath'));
    oldpath = path;
    cleanup = onCleanup(@() path(oldpath));
    addpath(fullfile(scriptDir,'..','..','help_funcs'));
    addpath(fullfile(scriptDir,'..','..','components'));
    addpath(fullfile(scriptDir,'..','..','..','verilog_support'));

    % Derive capacitor mismatch and initial thresholds using existing helpers
    caps = derrive_caps(N_bits, mismatch_LSC, LSC);
    [thresholds, ~] = cap_to_thr(caps, Vlow, Vhigh);

    % Generate a test input waveform in the ADC input range
    N = 2^N_bits - 1;
    f1 = (15.24532 / N) * fs;
    t = (0:cal_cycles-1)' / fs;
    amplitude = 0.8 * (Vhigh - Vlow) / 2;
    analog_in = ((Vhigh + Vlow)/2) + amplitude * sin(2*pi*f1*t);

    % Generate a matching dither sequence and apply 1 LSB amplitude
    dither_bits = uint8(gen_dither_sequence(cal_cycles));
    Vinc = (Vhigh - Vlow) / 2^N_bits;
    analog_in_dithered = analog_in + double(dither_bits) * Vinc;

    % Convert the dithered analog input to digital ADC codes
    data_in = uint16(flash_adc(analog_in_dithered, N_bits, Vhigh, Vlow, thresholds));

    % Provide a default DAC control input sequence for calibration
    DAC_ctrl_in = uint16(mod((0:cal_cycles-1)', 2^DAC_CTRL_WIDTH));

    if saveFile
        outPath = fullfile(scriptDir, 'calibration_input.mat');
        save(outPath, 'data_in', 'dither_bits', 'DAC_ctrl_in', 'thresholds', 'analog_in', 'Vhigh', 'Vlow', 'Vinc');

        % Write Verilog memory files for the testbench
        memDir = fullfile(scriptDir,'..','..','..','verilog','tb');
        if ~exist(memDir, 'dir')
            mkdir(memDir);
        end
        dataFile = fullfile(memDir, 'calibration_data_in.mem');
        ctrlFile = fullfile(memDir, 'calibration_DAC_ctrl_in.mem');

        dataWidthHex = ceil(N_bits / 4);
        ctrlWidthHex = ceil(DAC_CTRL_WIDTH / 4);

        fid = fopen(dataFile, 'w');
        fprintf(fid, ['%0' num2str(dataWidthHex) 'x\n'], data_in);
        fclose(fid);

        fid = fopen(ctrlFile, 'w');
        fprintf(fid, ['%0' num2str(ctrlWidthHex) 'x\n'], DAC_ctrl_in);
        fclose(fid);
    end
end

