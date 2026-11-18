%% Exercise 2: Frequency resolution, zero-padding and data length
% x[k] = cos(W0*k) + 0.6*cos(W1*k),  W0 = 12*pi/64,  W1 = 13*pi/64
%
% Key numbers:
%   Frequency gap        dW = W1 - W0 = pi/64
%   In DFT-bin units (L = 64) the two tones sit at bins 6 and 6.5.
%   Rectangular window of length N has main-lobe half width 2*pi/N.
%   To separate the two peaks we need roughly 2*pi/N <= dW  ->  N >= 128.
%   Zero padding (L > N) only samples the SAME DTFT more densely; it cannot
%   improve the physical resolution, which is set by the data length N.

clear; close all; clc;
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

W0 = 12*pi/64;
W1 = 13*pi/64;
sig = @(k) cos(W0*k) + 0.6*cos(W1*k);

%% ---------------- (1) N = 64, L = 64 / 128 / 512 ----------------
N  = 64;
k  = 0:N-1;
x  = sig(k);

Ls = [64 128 512];
Lref = 8192;                              % dense DTFT for reference
Xref = abs(fft(x, Lref));
wref = (0:Lref-1) * 2*pi / Lref;

fprintf('===== Exercise 2(1): N = 64 =====\n');
fprintf('Tone positions: W0/pi = %.4f, W1/pi = %.4f, gap = pi/64 = %.4f*pi\n', ...
        W0/pi, W1/pi, 1/64);
fprintf('Rectangular-window main-lobe half width 2*pi/N = %.4f*pi  (> gap -> not resolvable)\n', 2/N);

figure('Name', 'Ex2-1 N=64 zero padding', 'Position', [100 100 900 750]);
for i = 1:numel(Ls)
    L = Ls(i);
    X = abs(fft(x, L));
    w = (0:L-1) * 2*pi / L;
    subplot(numel(Ls), 1, i);
    plot(wref/pi, Xref, 'Color', [0.7 0.7 0.7]); hold on;
    stem(w/pi, X, 'filled', 'MarkerSize', 3);
    xline(W0/pi, 'r--'); xline(W1/pi, 'r--');
    xlim([0 0.5]); grid on;
    title(sprintf('N = 64, L = %d point FFT', L));
    xlabel('\Omega / \pi'); ylabel('|X[m]|');
    legend('DTFT (reference)', sprintf('%d-point DFT', L), '\Omega_0, \Omega_1', ...
           'Location', 'northeast');
    np = numel(localPeaks(X(1:floor(L/2)), 0.3*max(X)));
    fprintf('L = %4d : number of peaks found in [0, pi) = %d\n', L, np);
end
saveFig(gcf, figDir, 'Ex2_1_N64_zeropad');

%% ---------------- (2) Increase the data length N ----------------
% Choose N so that the two components are separated:
%   rectangular window : N >= 2*pi/dW = 128
%   Hamming window     : main-lobe half width 4*pi/N -> N >= 256
Ns = [64 128 256 512];
L  = 2048;                       % plenty of zero padding just for a smooth curve
fprintf('\n===== Exercise 2(2): increase the record length N (L = %d) =====\n', L);

figure('Name', 'Ex2-2 increase N (rectangular)', 'Position', [100 100 900 800]);
for i = 1:numel(Ns)
    N = Ns(i);
    x = sig(0:N-1);
    X = abs(fft(x, L));
    w = (0:L-1) * 2*pi / L;
    subplot(numel(Ns), 1, i);
    plot(w/pi, X, 'b', 'LineWidth', 1.1); hold on;
    xline(W0/pi, 'r--'); xline(W1/pi, 'r--');
    xlim([0.1 0.32]); grid on;
    title(sprintf('Rectangular window, N = %d, L = %d', N, L));
    xlabel('\Omega / \pi'); ylabel('|X(e^{j\Omega})|');
    pk = localPeaks(X(1:L/2), 0.3*max(X));
    fprintf('Rect, N = %4d : %d peak(s) at Omega/pi = %s\n', N, numel(pk), ...
            mat2str(round((pk-1)*2/L, 4)));
end
saveFig(gcf, figDir, 'Ex2_2_increaseN_rect');

%% ---------------- (2a) Best choice: N = L = 128 (no leakage) ----------------
% With N = 128, W0 = 2*pi*12/128 and W1 = 2*pi*13/128 fall EXACTLY on DFT bins
% 12 and 13 (an integer number of periods in the window), so a 128-point DFT
% without zero padding shows two clean lines and no leakage at all.
N = 128; Lb = 128;
x = sig(0:N-1);
X = abs(fft(x, Lb));
Xd = abs(fft(x, 4096));
figure('Name', 'Ex2-2a N=L=128', 'Position', [100 100 900 380]);
plot((0:4095)*2/4096, Xd, 'Color', [0.7 0.7 0.7]); hold on;
stem((0:Lb-1)*2/Lb, X, 'filled', 'MarkerSize', 4);
xline(W0/pi, 'r--'); xline(W1/pi, 'r--');
xlim([0 0.5]); grid on;
xlabel('\Omega / \pi'); ylabel('|X[m]|');
title('N = 128, L = 128: tones fall exactly on bins m = 12 and m = 13');
legend('DTFT of the 128-point record', '128-point DFT', '\Omega_0, \Omega_1');
[~, srt] = sort(X(1:Lb/2), 'descend');
fprintf('\nN = L = 128: two largest bins m = %s, |X| = %s (expected N/2 = 64 and 0.6*64 = 38.4)\n', ...
        mat2str(sort(srt(1:2)) - 1), mat2str(round(X(sort(srt(1:2))), 2)));
fprintf('All other bins in [0, pi): max |X| = %.2e  (no leakage)\n', ...
        max(X(setdiff(1:Lb/2, srt(1:2)))));
saveFig(gcf, figDir, 'Ex2_2a_N128_exactbins');

%% ---------------- (2b) Effect of window type ----------------
% Compare rectangular / Hann / Hamming / Blackman for several N.
% A tapered window lowers the side lobes (less leakage) but WIDENS the main
% lobe (Hann/Hamming: 4*pi/N, Blackman: 6*pi/N half width), so it needs a
% longer record to resolve two close lines.
winNames = {'Rectangular', 'Hann', 'Hamming', 'Blackman'};
Ncmp = [128 256 512];
fprintf('\n--- Window comparison (true positions 0.1875 and 0.2031) ---\n');
figure('Name', 'Ex2-3 window comparison', 'Position', [100 100 1000 900]);
for j = 1:numel(Ncmp)
    N = Ncmp(j);
    x = sig(0:N-1);
    subplot(numel(Ncmp), 1, j); hold on;
    for i = 1:numel(winNames)
        win = makeWindow(winNames{i}, N);
        X = abs(fft(x .* win, L));
        X = X / max(X);
        w = (0:L-1) * 2*pi / L;
        plot(w/pi, 20*log10(X + eps), 'LineWidth', 1.1);
        pk = localPeaks(X(1:L/2), 0.3);
        fprintf('%-11s N = %3d : %d peak(s) at Omega/pi = %s\n', winNames{i}, N, ...
                numel(pk), mat2str(round((pk-1)*2/L, 4)));
    end
    xline(W0/pi, 'k--'); xline(W1/pi, 'k--');
    xlim([0.1 0.32]); ylim([-80 5]); grid on;
    title(sprintf('Window comparison, N = %d, L = %d (normalized, dB)', N, L));
    xlabel('\Omega / \pi'); ylabel('Magnitude (dB)');
    legend([winNames, {'\Omega_0, \Omega_1'}], 'Location', 'southwest');
end
saveFig(gcf, figDir, 'Ex2_3_windows');

fprintf(['\nConclusion: with N = 64 the two lines cannot be separated whatever L is,\n' ...
         'because zero-padding does not add information. N = L = 128 (rectangular,\n' ...
         'integer periods) resolves them exactly; tapered windows need N >= ~384.\n' ...
         'Note: at N = 128 the tapered windows show two peaks at WRONG positions,\n' ...
         'caused by interference of the two overlapping main lobes.\n']);

%% ======================= local functions =======================
function idx = localPeaks(y, thr)
% Simple local-maximum peak picker (no toolbox needed)
y = y(:).';
idx = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end) & y(2:end-1) > thr) + 1;
end

function w = makeWindow(name, N)
n = 0:N-1;
switch lower(name)
    case 'rectangular', w = ones(1, N);
    case 'hann',        w = 0.5 - 0.5*cos(2*pi*n/(N-1));
    case 'hamming',     w = 0.54 - 0.46*cos(2*pi*n/(N-1));
    case 'blackman',    w = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1));
end
end

function saveFig(fig, figDir, name)
exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 200);
end
