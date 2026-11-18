%% Exercise 4: Spectral components of a chord (hexian1.wav)
% (1) Use the DFT to find which notes make up the chord.
% (2) Whole note = 0.2 s  ->  16th note = 0.2/16 = 0.0125 s.
%     Theory: a record of length T gives a frequency resolution of about 1/T
%     (rectangular window). 1/0.0125 s = 80 Hz, but neighbouring notes of the
%     chord are only ~70-76 Hz apart, so their main lobes overlap and the DFT
%     cannot separate them, whatever the zero padding. Verified below, and a
%     model-based remedy (least squares on the note grid) is demonstrated.

clear; close all; clc;
here   = fileparts(mfilename('fullpath'));
figDir = fullfile(here, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
noteNames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};

%% ======================= (1) chord analysis =======================
[x, fs] = audioread(fullfile(here, 'hexian1.wav'));
x = x(:, 1);
N = length(x);  T = N/fs;  t = (0:N-1)/fs;
L = 4 * 2^nextpow2(N);                       % zero padding for a smooth plot

fprintf('===== Exercise 4(1): hexian1.wav =====\n');
fprintf('fs = %d Hz, N = %d, T = %.3f s, resolution fs/N = %.2f Hz, L = %d\n', ...
        fs, N, T, fs/N, L);

X = abs(fft(x, L)) / N * 2;
f = (0:L-1) * fs / L;
half = 1:L/2;
[pkVal, pkLoc] = findpeaks(X(half), 'MinPeakHeight', 0.3*max(X), ...
                           'MinPeakDistance', round(10/(fs/L)));
fChord = f(pkLoc);
fprintf('Spectral peaks of the chord:\n');
chordNames = cell(1, numel(fChord));
for i = 1:numel(fChord)
    [chordNames{i}, ~, c] = freq2note(fChord(i));
    fprintf('  %7.2f Hz  amp %.3f  -> %-4s (%+5.1f cents)\n', fChord(i), pkVal(i), chordNames{i}, c);
end
fprintf('==> The chord is %s  (D major triad: root D, major third F#, perfect fifth A)\n', ...
        strjoin(chordNames, ' + '));

figure('Name', 'Ex4-1 chord spectrum', 'Position', [100 100 1000 720]);
subplot(3,1,1);
plot(t, x); grid on; xlim([0 T]);
xlabel('t (s)'); ylabel('x(t)'); title(sprintf('hexian1.wav, f_s = %d Hz, N = %d, T = %.1f s', fs, N, T));
subplot(3,1,2);
plot(t*1000, x); grid on; xlim([0 30]);
xlabel('t (ms)'); ylabel('x(t)'); title('First 30 ms of the waveform');
subplot(3,1,3);
plot(f(half), X(half)); hold on; grid on;
plot(fChord, pkVal, 'rv', 'MarkerFaceColor', 'r');
for i = 1:numel(fChord)
    text(fChord(i), pkVal(i)*1.08, sprintf('%.1f Hz\n%s', fChord(i), chordNames{i}), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
xlim([200 550]); ylim([0 max(pkVal)*1.45]);
xlabel('f (Hz)'); ylabel('|X|'); title(sprintf('%d-point DFT of the chord', L));
saveFig(gcf, figDir, 'Ex4_1_chord_spectrum');

%% ======================= (2) 16th-note problem =======================
Twhole = 0.2;
durs   = Twhole ./ [16 8 4 2 1];             % 1/16, 1/8, 1/4, 1/2, whole note
labels = {'1/16 note','1/8 note','1/4 note','1/2 note','whole note'};
dfMin  = min(diff(sort(fChord)));
fprintf('\n===== Exercise 4(2): short notes =====\n');
fprintf('16th note duration = %.4f s -> %d samples at fs = %d Hz\n', durs(1), round(durs(1)*fs), fs);
fprintf('Smallest spacing between chord notes = %.1f Hz\n', dfMin);
fprintf('Rule of thumb: need T >= 1/df = %.1f ms (rectangular), >= 2/df = %.1f ms (Hann)\n', ...
        1000/dfMin, 2000/dfMin);

Lz = 8192;
fz = (0:Lz-1) * fs / Lz;
figure('Name', 'Ex4-2 duration vs resolution', 'Position', [100 100 1000 950]);
fprintf('\n  duration    N   fs/N(Hz) | rect: peaks (Hz) -> notes            | Hann: peaks (Hz) -> notes\n');
for i = 1:numel(durs)
    Ns = round(durs(i) * fs);
    xs = x(1:Ns);
    Xr = abs(fft(xs, Lz));            Xr = Xr / max(Xr);
    Xh = abs(fft(xs .* hann(Ns), Lz)); Xh = Xh / max(Xh);
    [~, locR] = findpeaks(Xr(1:Lz/2), 'MinPeakHeight', 0.5);
    [~, locH] = findpeaks(Xh(1:Lz/2), 'MinPeakHeight', 0.5);
    nmR = strjoin(arrayfun(@(v) freq2note(v), fz(locR), 'UniformOutput', false), ' ');
    nmH = strjoin(arrayfun(@(v) freq2note(v), fz(locH), 'UniformOutput', false), ' ');
    fprintf('  %-10s %4d %7.1f  | %-20s -> %-13s | %-20s -> %s\n', labels{i}, Ns, fs/Ns, ...
            mat2str(round(fz(locR))), nmR, mat2str(round(fz(locH))), nmH);
    subplot(numel(durs), 1, i);
    plot(fz(1:Lz/2), Xr(1:Lz/2), 'b', 'LineWidth', 1.1); hold on;
    plot(fz(1:Lz/2), Xh(1:Lz/2), 'Color', [0.85 0.4 0.1]);
    for k = 1:numel(fChord), xline(fChord(k), 'k--'); end
    xlim([100 650]); ylim([0 1.1]); grid on;
    title(sprintf('%s: T = %.4f s, N = %d, f_s/N = %.1f Hz', labels{i}, durs(i), Ns, fs/Ns));
    ylabel('|X| (norm.)');
    if i == 1, legend('Rectangular', 'Hann', 'true notes', 'Location', 'northeast'); end
end
xlabel('f (Hz)');
saveFig(gcf, figDir, 'Ex4_2_duration_vs_resolution');

% Zero padding does not help: the 16th-note record with L = 128 ... 16384
figure('Name', 'Ex4-3 16th note zero padding', 'Position', [100 100 1000 600]);
Ns = round(durs(1) * fs);
xs = x(1:Ns);
subplot(2,1,1);
plot((0:Ns-1)/fs*1000, xs, '.-'); grid on;
xlabel('t (ms)'); ylabel('x'); title(sprintf('A 16th note of the chord: %d samples (%.1f ms)', Ns, durs(1)*1000));
subplot(2,1,2); hold on; grid on;
Lp = [128 1024 16384];
for i = 1:numel(Lp)
    Xp = abs(fft(xs, Lp(i)));
    plot((0:Lp(i)/2-1)*fs/Lp(i), Xp(1:Lp(i)/2), '.-', 'MarkerSize', 8 - 2*i);
end
for k = 1:numel(fChord), xline(fChord(k), 'k--'); end
xlim([0 1000]);
xlabel('f (Hz)'); ylabel('|X[m]|');
legend('L = 128', 'L = 1024', 'L = 16384', 'true notes');
title('16th note: zero padding only smooths the curve - the peaks stay at WRONG frequencies');
saveFig(gcf, figDir, 'Ex4_3_16th_zero_padding');

% Same experiment on a synthesized C major chord (C4 E4 G4) to show that
% the problem is general and not specific to this recording.
fC = [261.63 329.63 392];
tt = (0:round(durs(1)*fs)-1)/fs;
xc = sum(cos(2*pi*fC(:)*tt), 1).';
tl = (0:fs-1)/fs;
xcl = sum(cos(2*pi*fC(:)*tl), 1).';
figure('Name', 'Ex4-4 synthetic C chord', 'Position', [100 100 1000 420]);
Xc  = abs(fft(xc, Lz));  Xc  = Xc/max(Xc);
Xcl = abs(fft(xcl, Lz)); Xcl = Xcl/max(Xcl);
plot(fz(1:Lz/2), Xcl(1:Lz/2), 'Color', [0.6 0.6 0.6]); hold on;
plot(fz(1:Lz/2), Xc(1:Lz/2), 'b', 'LineWidth', 1.3);
for k = 1:3, xline(fC(k), 'r--'); end
xlim([100 600]); grid on;
xlabel('f (Hz)'); ylabel('|X| (norm.)');
legend('T = 1 s', 'T = 12.5 ms (16th note)', 'C4, E4, G4');
title('Synthetic C-major chord: the 16th note cannot be resolved by the DFT');
saveFig(gcf, figDir, 'Ex4_4_synthetic_C_chord');

%% ---------- Remedy: model-based least-squares fit on the note grid ----------
% Music notes lie on a known grid f = 440*2^((m-69)/12). For every set of 3
% candidate notes we fit  x[n] = sum_k a_k cos(2*pi*f_k*n/fs) + b_k sin(...)
% by least squares and keep the set with the smallest residual. This uses
% prior knowledge instead of DFT resolution, so 100 samples are enough.
cand = 60:83;                                 % C4 ... B5
fprintf('\n--- Least-squares note-grid fit on ONE 16th note (%d samples) ---\n', Ns);
testSigs = {xs, 'hexian1 16th note (clean)'; ...
            xs + 10^(-20/20)*rms(xs)*randn(size(xs)), 'hexian1 16th note + noise (SNR 20 dB)'; ...
            xc, 'synthetic C4 E4 G4 16th note'};
rng(0);
for s = 1:size(testSigs, 1)
    [best, res] = lsChordFit(testSigs{s,1}, fs, cand, 3);
    nm = arrayfun(@(m) sprintf('%s%d', noteNames{mod(m,12)+1}, floor(m/12)-1), best, 'UniformOutput', false);
    fprintf('  %-40s -> %s   (relative residual %.2e)\n', testSigs{s,2}, strjoin(nm, ' '), res);
end

%% ======================= local functions =======================
function [bestSet, bestRes] = lsChordFit(x, fs, cand, K)
x = x(:);
n = (0:numel(x)-1).';
combos = nchoosek(cand, K);
bestRes = inf; bestSet = [];
for i = 1:size(combos, 1)
    fk = 440 * 2.^((combos(i,:) - 69)/12);
    A = [cos(2*pi*n*fk/fs), sin(2*pi*n*fk/fs)];
    c = A \ x;
    r = norm(x - A*c)^2 / norm(x)^2;
    if r < bestRes, bestRes = r; bestSet = combos(i,:); end
end
end

function [name, midi, cents] = freq2note(f)
names = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
m = 69 + 12*log2(f/440);
midi = round(m);
cents = 100*(m - midi);
name = sprintf('%s%d', names{mod(midi,12)+1}, floor(midi/12) - 1);
end

function saveFig(fig, figDir, name)
exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 200);
end
