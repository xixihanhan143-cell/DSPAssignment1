%% Exercise 5: Spectral analysis of a scale with harmonics (yueyin2.wav)

clear; close all; clc;
here   = fileparts(mfilename('fullpath'));
figDir = fullfile(here, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

[x, fs] = audioread(fullfile(here, 'yueyin2.wav'));
x = x(:, 1);
N = length(x);  T = N/fs;  t = (0:N-1)/fs;
fprintf('===== Exercise 5: yueyin2.wav =====\n');
fprintf('fs = %d Hz, N = %d, T = %.2f s\n', fs, N, T);

%% Part A: whole-record DFT
L = 2^nextpow2(N);
X = abs(fft(x, L)) / N * 2;
f = (0:L-1) * fs / L;
half = 1:L/2;
[pkVal, pkLoc] = findpeaks(X(half), 'MinPeakHeight', 0.05*max(X), ...
                           'MinPeakDistance', round(10/(fs/L)));
fPk = f(pkLoc).';
fprintf('\nPart A: the whole-record DFT shows %d spectral lines (0 - %d Hz):\n', numel(fPk), fs/2);

scaleF = [261.63 293.66 329.63 349.23 392 440 493.88 523.25];
scaleN = {'C4','D4','E4','F4','G4','A4','B4','C5'};
fprintf('   f (Hz)   rel.amp   possible sources within 3 Hz (note x harmonic number)\n');
nAmb = 0;
for i = 1:numel(fPk)
    src = {};
    for k = 1:numel(scaleF)
        n = round(fPk(i)/scaleF(k));
        if n >= 1 && abs(n*scaleF(k) - fPk(i)) < 3
            src{end+1} = sprintf('%s x%d', scaleN{k}, n); %#ok<SAGROW>
        end
    end
    if numel(src) > 1, nAmb = nAmb + 1; flag = '  <-- ambiguous'; else, flag = ''; end
    fprintf('  %7.1f   %6.3f    %s%s\n', fPk(i), pkVal(i)/max(pkVal), strjoin(src, ', '), flag);
end
fprintf('%d of %d lines can come from more than one note.\n', nAmb, numel(fPk));

figure('Name', 'Ex5-1 whole spectrum', 'Position', [100 100 1000 650]);
subplot(2,1,1);
plot(t, x); grid on; xlim([0 T]);
xlabel('t (s)'); ylabel('x(t)'); title('yueyin2.wav (scale with harmonics)');
subplot(2,1,2);
plot(f(half), X(half)); hold on; grid on;
plot(fPk, pkVal, 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
for i = 1:numel(fPk)
    text(fPk(i), pkVal(i) + 0.004, sprintf('%.0f', fPk(i)), 'Rotation', 90, 'FontSize', 7);
end
xlim([0 2300]); ylim([0 max(pkVal)*1.35]);
xlabel('f (Hz)'); ylabel('|X|');
title(sprintf('%d-point DFT of the whole record: %d lines, notes cannot be told apart', L, numel(fPk)));
saveFig(gcf, figDir, 'Ex5_1_whole_spectrum');

x1 = audioread(fullfile(here, 'yueyin1.wav'));
X1 = abs(fft(x1(:,1), L)) / N * 2;
figure('Name', 'Ex5-2 compare yueyin1 yueyin2', 'Position', [100 100 1000 420]);
plot(f(half), X1(half), 'Color', [0.6 0.6 0.6]); hold on;
plot(f(half), X(half), 'b');
grid on; xlim([0 2300]);
xlabel('f (Hz)'); ylabel('|X|');
legend('yueyin1 (pure tones)', 'yueyin2 (with harmonics)');
title('Whole-record spectra of the two scale signals');
saveFig(gcf, figDir, 'Ex5_2_compare_yueyin1_2');

%% Part B: per-note analysis
figure('Name', 'Ex5-3 spectrogram', 'Position', [100 100 1000 450]);
spectrogram(x, hann(1024), 768, 4096, fs, 'yaxis');
ylim([0 2.3]); colormap(jet); clim([-90 -10]);
xlabel('Time (s)'); ylabel('Frequency (kHz)');
cb = colorbar; cb.Label.String = 'Power/frequency (dB/Hz)';
title('Spectrogram of yueyin2.wav: each note = a ladder of harmonics n f_0');
saveFig(gcf, figDir, 'Ex5_3_spectrogram');

res = analyseHarmonics(x, fs);
nSeg = numel(res);
printTable(res, fs);

figure('Name', 'Ex5-4 per-note spectra', 'Position', [100 100 1150 800]);
for k = 1:nSeg
    subplot(4, 2, k);
    plot(res(k).f, res(k).X); hold on; grid on;
    hv = find(res(k).A > 0.01*res(k).A(1));          % harmonics actually present
    plot(hv*res(k).f0, res(k).A(hv), 'ro', 'MarkerSize', 5);
    for n = hv
        text(n*res(k).f0, res(k).A(n) + 0.05, sprintf('%d', n), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    xlim([0 2300]); ylim([0 0.65]);
    title(sprintf('Note %d: f_0 = %.2f Hz (%s)', k, res(k).f0, res(k).name));
    if k > nSeg - 2, xlabel('f (Hz)'); end
    ylabel('A');
end
saveFig(gcf, figDir, 'Ex5_4_per_note_spectra');

Mmax = max(arrayfun(@(r) numel(r.A), res));
Amat = zeros(nSeg, Mmax);
for k = 1:nSeg, Amat(k, 1:numel(res(k).A)) = res(k).A / res(k).A(1); end
Mshow = find(any(Amat > 0.005, 1), 1, 'last');
figure('Name', 'Ex5-5 harmonic amplitudes', 'Position', [100 100 900 400]);
bar(Amat(:, 1:Mshow + 1)); grid on;
set(gca, 'XTickLabel', {res.name});
legend(arrayfun(@(n) sprintf('A_%d/A_1', n), 1:Mshow + 1, 'UniformOutput', false), 'Location', 'northeastoutside');
ylabel('relative amplitude'); title('Estimated harmonic composition of every note');
saveFig(gcf, figDir, 'Ex5_5_harmonic_amplitudes');

%% Part C: verification
fprintf('\nC1: whole-record line heights predicted from the per-note results\n');
fprintf('   f (Hz)   measured   predicted   contributions\n');
ref = median(pkVal(arrayfun(@(i) numel(contrib(fPk(i), res)) == 1 && ...
         contrib(fPk(i), res).n == 1, 1:numel(fPk)))) / mean(arrayfun(@(r) r.A(1), res));
errs = zeros(numel(fPk), 1);
for i = 1:numel(fPk)
    c = contrib(fPk(i), res);
    pred = sum([c.A]) * ref;
    errs(i) = (pkVal(i) - pred) / pred;
    fprintf('  %7.1f   %7.4f    %7.4f    %s\n', fPk(i), pkVal(i), pred, ...
            strjoin(arrayfun(@(s) sprintf('%s x%d', s.name, s.n), c, 'UniformOutput', false), ' + '));
end
fprintf('Mean |relative error| over all lines: %.1f %%\n', 100*mean(abs(errs)));
fprintf(['(Shared lines are slightly LOWER than the plain sum: the contributions come from\n' ...
         ' different time intervals and add as complex numbers with different phases.)\n']);

xr = zeros(N, 1);
for k = 1:nSeg
    idx = res(k).seg(1):res(k).seg(2);
    tk  = (0:numel(idx)-1).' / fs;
    w   = hann(numel(idx));
    env = noteEnvelope(x(idx), fs, res(k).f0);
    env = env / (sum(w .* env) / sum(w));
    M   = numel(res(k).A);
    B   = [cos(2*pi*tk*(1:M)*res(k).f0), sin(2*pi*tk*(1:M)*res(k).f0)] .* env;
    cLS = B \ x(idx);
    phi = atan2(-cLS(M+1:end), cLS(1:M)).';
    tone = cos(2*pi*tk*(1:M)*res(k).f0 + phi) * res(k).A.';
    xr(idx) = tone .* env;
end
Xr = abs(fft(xr, L)) / N * 2;
c = corrcoef(X(half), Xr(half));
snr = 10*log10(sum(x.^2) / sum((x - xr).^2));
fprintf('\nC2: resynthesis - spectral correlation = %.4f, waveform SNR = %.1f dB\n', c(1,2), snr);
audiowrite(fullfile(here, 'yueyin2_resynth.wav'), xr / max(abs(xr)) * 0.9, fs);
fprintf('    resynthesised audio written to yueyin2_resynth.wav\n');

figure('Name', 'Ex5-6 resynthesis', 'Position', [100 100 1000 650]);
subplot(2,1,1);
plot(t, x, 'Color', [0.6 0.6 0.6]); hold on; plot(t, xr, 'b:');
xlim([0.95 1.05]); grid on;
legend('original', 'resynthesised'); xlabel('t (s)'); ylabel('x');
title(sprintf('Waveform around a note boundary (SNR = %.1f dB)', snr));
subplot(2,1,2);
plot(f(half), X(half), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); hold on;
plot(f(half), Xr(half), 'b:');
xlim([0 2300]); grid on;
legend('original', 'resynthesised'); xlabel('f (Hz)'); ylabel('|X|');
title(sprintf('Whole-record spectra (correlation %.4f)', c(1,2)));
saveFig(gcf, figDir, 'Ex5_6_resynthesis');

rng(1);
gMidi  = [67 69 71 72 74 76 78 79];                 % G4 A4 B4 C5 D5 E5 F#5 G5
f0True = 440 * 2.^((gMidi - 69)/12);
Atrue  = [1 0.6 0.35 0.15 0.05];
durs   = [0.4 0.4 0.6 0.4 0.4 0.6 0.4 0.8];
xs = [];
for k = 1:numel(gMidi)
    nk = round(durs(k)*fs);  tk = (0:nk-1).'/fs;
    tone = zeros(nk, 1);
    for n = 1:numel(Atrue), tone = tone + Atrue(n)*cos(2*pi*n*f0True(k)*tk + 2*pi*rand); end
    xs = [xs; 0.5 * tone .* adsr(nk, fs)]; %#ok<AGROW>
end
xs = xs + 10^(-30/20) * rms(xs) * randn(size(xs));
rs = analyseHarmonics(xs, fs);
fprintf('\nC3: synthetic G-major scale with known harmonics A = %s\n', mat2str(Atrue));
fprintf('  note  true f0   est f0   err(cent)   estimated A_n/A_1\n');
maxAerr = 0; maxCent = 0;
for k = 1:min(numel(rs), numel(f0True))
    Ae = rs(k).A(1:numel(Atrue)) / rs(k).A(1);
    cent = 1200*log2(rs(k).f0/f0True(k));
    maxAerr = max(maxAerr, max(abs(Ae - Atrue)));
    maxCent = max(maxCent, abs(cent));
    fprintf('  %-4s %8.2f %8.2f   %+6.2f     %s\n', rs(k).name, f0True(k), rs(k).f0, cent, mat2str(round(Ae, 3)));
end
fprintf('  detected %d of %d notes, max pitch error %.2f cent, max A_n error %.3f\n', ...
        numel(rs), numel(f0True), maxCent, maxAerr);

figure('Name', 'Ex5-7 synthetic test', 'Position', [100 100 900 380]);
Ae = cell2mat(arrayfun(@(r) r.A(1:numel(Atrue)) / r.A(1), rs(:), 'UniformOutput', false));
bar(Ae.'); hold on;
plot(1:numel(Atrue), Atrue, 'k*-', 'LineWidth', 1.5, 'MarkerSize', 9);
grid on; ylim([0 1.12]); xlabel('harmonic number n'); ylabel('A_n / A_1');
legend([{rs.name}, {'true values'}], 'Location', 'northeastoutside');
title('Synthetic test: estimated harmonic amplitudes of all 8 notes vs ground truth');
saveFig(gcf, figDir, 'Ex5_7_synthetic_test');

%% local functions
function res = analyseHarmonics(x, fs)
seg = segmentNotes(x, fs);
Lseg = 16384;
fAx = (0:Lseg/2-1).' * fs / Lseg;
res = struct('seg', {}, 'f0', {}, 'name', {}, 'A', {}, 'phi', {}, 'f', {}, 'X', {});
for k = 1:size(seg, 1)
    xs = x(seg(k,1):seg(k,2));
    w = hann(numel(xs));
    Xc = fft(xs .* w, Lseg);
    Xm = abs(Xc(1:Lseg/2)) * 2 / sum(w);              % amplitude-calibrated
    [pv, pl] = findpeaks(Xm, 'MinPeakHeight', 0.25*max(Xm), ...
                         'MinPeakDistance', round(20/(fs/Lseg)));
    pl = pl(fAx(pl) > 60);
    f0 = refinePeak(Xm, pl(1)) * fs / Lseg;
    M = floor((fs/2 - 50) / f0);
    A = zeros(1, M); phi = zeros(1, M);
    for n = 1:M
        band = find(abs(fAx - n*f0) <= 0.02*f0);
        [A(n), j] = max(Xm(band));
        phi(n) = angle(Xc(band(j)));
    end
    res(k).seg = seg(k,:); res(k).f0 = f0;  res(k).name = freq2note(f0);
    res(k).A = A;  res(k).phi = phi;  res(k).f = fAx;  res(k).X = Xm;
end
end

function printTable(res, fs)
fprintf('\nPart B: per-note analysis (segment -> Hann-windowed DFT -> f0 -> A_n)\n');
fprintf('  note  time (s)       f0 (Hz)   name   A_1     A_2/A_1  A_3/A_1  A_4/A_1  A_5/A_1  harmonics\n');
for k = 1:numel(res)
    r = res(k);
    A = r.A / r.A(1);
    sig = find(A > 0.01);
    fprintf('  %2d   %.3f-%.3f  %8.2f   %-5s  %.3f   %6.3f   %6.3f   %6.3f   %6.3f   %s Hz\n', ...
            k, (r.seg(1)-1)/fs, r.seg(2)/fs, r.f0, r.name, r.A(1), A(2), A(3), A(4), A(5), ...
            mat2str(round(sig * r.f0, 1)));
end
end

function c = contrib(fp, res)
c = struct('name', {}, 'n', {}, 'A', {});
for k = 1:numel(res)
    n = round(fp / res(k).f0);
    if n >= 1 && n <= numel(res(k).A) && abs(n*res(k).f0 - fp) < 3 && res(k).A(n) > 0.01*res(k).A(1)
        c(end+1) = struct('name', res(k).name, 'n', n, 'A', res(k).A(n)); %#ok<AGROW>
    end
end
end

function env = noteEnvelope(xs, fs, f0)
win = max(8, round(2*fs/f0));
env = sqrt(movmean(xs.^2, win)) * sqrt(2);
end

function e = adsr(n, fs)
a = round(0.05*fs); d = round(0.05*fs); r = round(0.1*fs);
s = max(n - a - d - r, 0);
e = [linspace(0,1,a), linspace(1,0.7,d), 0.7*ones(1,s), linspace(0.7,0,r)].';
e = e(1:n);
end

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
note(en < 0.01*max(en)) = NaN;
seg = [];
i = 1;
while i <= nFr
    if isnan(note(i)), i = i + 1; continue; end
    j = i;
    while j < nFr && note(j+1) == note(i), j = j + 1; end
    if (j - i + 1) * hop / fs > 0.05
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

function p = refinePeak(X, k)
if k <= 1 || k >= numel(X), p = k - 1; return; end
a = X(k-1); b = X(k); c = X(k+1);
p = (k - 1) + 0.5*(a - c)/(a - 2*b + c);
end

function name = freq2note(f)
names = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
m = round(69 + 12*log2(f/440));
name = sprintf('%s%d', names{mod(m,12)+1}, floor(m/12) - 1);
end

function saveFig(fig, figDir, name)
exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 200);
end
