%% Exercise 3: Spectrum of the scale signal yueyin1.wav and its key

clear; close all; clc;
here   = fileparts(mfilename('fullpath'));
figDir = fullfile(here, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% 1) basic parameters
[x, fs] = audioread(fullfile(here, 'yueyin1.wav'));
x = x(:, 1);                      % mono
N = length(x);
T = N / fs;
t = (0:N-1) / fs;

fprintf('===== Exercise 3: yueyin1.wav =====\n');
fprintf('Sampling frequency fs = %d Hz\n', fs);
fprintf('Number of samples  N  = %d\n', N);
fprintf('Duration           T  = N/fs = %.3f s\n', T);

L = 2^nextpow2(N);
fprintf('DFT length L = %d (>= N), bin spacing fs/L = %.4f Hz, physical resolution fs/N = %.4f Hz\n', ...
        L, fs/L, fs/N);

%% 2) DFT of the whole signal
X  = abs(fft(x, L)) / N * 2;      % single-sided amplitude
f  = (0:L-1) * fs / L;
half = 1:L/2;

[pkVal, pkLoc] = findpeaks(X(half), 'MinPeakHeight', 0.2*max(X), ...
                           'MinPeakDistance', round(10/(fs/L)));
fPk = f(pkLoc);
fprintf('\nPeaks in the spectrum of the whole signal:\n');
for i = 1:numel(fPk)
    [nm, ~, cents] = freq2note(fPk(i));
    fprintf('  %8.2f Hz   amp %.3f   -> %-4s (%+5.1f cents)\n', fPk(i), pkVal(i), nm, cents);
end

refF = [261.63 293.66 329.63 349.23 392 440 493.88 523.25];
figure('Name', 'Ex3-1 whole spectrum', 'Position', [100 100 1000 700]);
subplot(3,1,1);
plot(t, x); grid on; xlim([0 T]);
xlabel('t (s)'); ylabel('x(t)'); title(sprintf('yueyin1.wav, f_s = %d Hz, N = %d, T = %.1f s', fs, N, T));
subplot(3,1,2);
plot(f(half), X(half)); grid on; xlim([0 fs/2]);
xlabel('f (Hz)'); ylabel('|X|'); title(sprintf('%d-point DFT of the whole signal (0 \\sim f_s/2)', L));
subplot(3,1,3);
plot(f(half), X(half)); hold on; grid on;
plot(fPk, pkVal, 'rv', 'MarkerFaceColor', 'r');
for i = 1:numel(fPk)
    text(fPk(i), pkVal(i)*1.08, sprintf('%.1f', fPk(i)), 'HorizontalAlignment', 'center', 'FontSize', 8);
end
xlim([230 560]); ylim([0 max(pkVal)*1.25]);
xlabel('f (Hz)'); ylabel('|X|'); title('Zoom 230 - 560 Hz: eight spectral lines');
saveFig(gcf, figDir, 'Ex3_1_whole_spectrum');

%% 3) short-time analysis: order of the notes
win = 1024; hop = 256; nfft = 4096;
figure('Name', 'Ex3-2 spectrogram', 'Position', [100 100 1000 420]);
spectrogram(x, hann(win), win-hop, nfft, fs, 'yaxis');
ylim([0 1]); colormap(jet); clim([-80 -10]);
xlabel('Time (s)'); ylabel('Frequency (kHz)');
cb = colorbar; cb.Label.String = 'Power/frequency (dB/Hz)';
title('Spectrogram of yueyin1.wav (Hann 1024, hop 256)');
saveFig(gcf, figDir, 'Ex3_2_spectrogram');

seg = segmentNotes(x, fs);
nSeg = size(seg, 1);
fprintf('\nDetected %d notes (pitch-change segmentation):\n', nSeg);

%% 4) DFT of every note
Lseg = 8192;
f0 = zeros(nSeg, 1); names = cell(nSeg, 1); midi = zeros(nSeg, 1);
figure('Name', 'Ex3-3 per-note spectra', 'Position', [100 100 1100 750]);
for k = 1:nSeg
    xs = x(seg(k,1):seg(k,2));
    Ns = numel(xs);
    Xs = abs(fft(xs .* hann(Ns), Lseg));
    fs_ = (0:Lseg-1) * fs / Lseg;
    [~, im] = max(Xs(1:Lseg/2));
    f0(k) = parabolicPeak(Xs, im) * fs / Lseg;     % refine between bins
    [names{k}, midi(k), cents] = freq2note(f0(k));
    fprintf('  note %d: %.3f-%.3f s (%4d samples), f = %7.2f Hz -> %-4s (%+5.1f cents)\n', ...
            k, (seg(k,1)-1)/fs, seg(k,2)/fs, Ns, f0(k), names{k}, cents);
    subplot(2, ceil(nSeg/2), k);
    plot(fs_(1:Lseg/2), Xs(1:Lseg/2) / max(Xs)); grid on;
    xlim([200 600]); ylim([0 1.1]);
    xline(f0(k), 'r--');
    title(sprintf('Note %d: %.1f Hz (%s)', k, f0(k), names{k}));
    xlabel('f (Hz)');
end
saveFig(gcf, figDir, 'Ex3_3_per_note_spectra');

%% 5) key determination
steps = diff(midi).';
fprintf('\nInterval pattern (semitones): %s\n', mat2str(steps));
major = [2 2 1 2 2 2 1];
if isequal(steps, major)
    tonic = regexprep(names{1}, '\d', '');
    fprintf('Pattern = W W H W W W H  -> MAJOR scale, tonic (do, "1") = %s (%.2f Hz)\n', ...
            names{1}, f0(1));
    fprintf('==> The scale is in %s major (1 = %s).\n', tonic, tonic);
else
    fprintf('Pattern is not a plain major scale - check manually.\n');
end

fprintf('\n  Note   measured(Hz)   table(Hz)   error(Hz)\n');
solfa = {'1','2','3','4','5','6','7','i'};
for k = 1:min(nSeg, numel(refF))
    fprintf('  %-4s   %10.2f   %9.2f   %+8.2f\n', solfa{k}, f0(k), refF(k), f0(k)-refF(k));
end

figure('Name', 'Ex3-4 waveform with notes', 'Position', [100 100 1000 330]);
plot(t, x, 'Color', [0.3 0.3 0.8]); hold on; grid on; xlim([0 T]);
yl = ylim;
for k = 1:nSeg
    xline((seg(k,1)-1)/fs, 'k:');
    text(mean(seg(k,:))/fs, yl(2)*0.9, sprintf('%s\n%.1f Hz', names{k}, f0(k)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end
xlabel('t (s)'); ylabel('x(t)'); title('yueyin1.wav segmented into notes');
saveFig(gcf, figDir, 'Ex3_4_waveform_notes');

%% local functions
function seg = segmentNotes(x, fs)
win = 512; hop = 128; nfft = 4096;
nFr = floor((numel(x) - win) / hop) + 1;
w = hann(win);
note = zeros(nFr, 1); en = zeros(nFr, 1);
fAx = (0:nfft/2-1) * fs / nfft;
band = fAx > 80 & fAx < 2000;
for i = 1:nFr
    fr = x((i-1)*hop + (1:win)) .* w;
    F = abs(fft(fr, nfft)); F = F(1:nfft/2);
    en(i) = sum(fr.^2);
    Fb = F; Fb(~band) = 0;
    [~, im] = max(Fb);
    note(i) = round(69 + 12*log2(fAx(im)/440));
end
note(en < 0.01*max(en)) = NaN;              % silence
seg = [];
i = 1;
while i <= nFr
    if isnan(note(i)), i = i + 1; continue; end
    j = i;
    while j < nFr && note(j+1) == note(i), j = j + 1; end
    if (j - i + 1) * hop / fs > 0.05            % ignore very short runs
        seg(end+1, :) = [i j]; %#ok<AGROW>
    end
    i = j + 1;
end
sStart = (seg(:,1) - 1) * hop + win/2;
sStart(1) = 1;
seg = [sStart, [sStart(2:end) - 1; numel(x)]];
e = movmean(x.^2, round(0.005*fs));
r = round(0.04*fs);
for b = 2:size(seg, 1)
    lo = max(1, seg(b,1) - r); hi = min(numel(x), seg(b,1) + r);
    [~, m] = min(e(lo:hi));
    seg(b,1) = lo + m - 1;  seg(b-1,2) = seg(b,1) - 1;
end
end

function [name, midi, cents] = freq2note(f)
names = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
m = 69 + 12*log2(f/440);
midi = round(m);
cents = 100*(m - midi);
name = sprintf('%s%d', names{mod(midi,12)+1}, floor(midi/12) - 1);
end

function p = parabolicPeak(X, k)
if k <= 1 || k >= numel(X), p = k - 1; return; end
a = X(k-1); b = X(k); c = X(k+1);
p = (k - 1) + 0.5*(a - c)/(a - 2*b + c);
end

function saveFig(fig, figDir, name)
exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 200);
end
