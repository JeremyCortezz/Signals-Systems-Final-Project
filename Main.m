clear; clc; close all;
%% ---- Settings -----------------------------------------------------------
filename  = ['S003_PSG_df.csv'];
% Detection
winSec =  1;
maxLagSec = 40;
lagOverrideSec = [];             % [] = estimate; or force e.g. 8
thrEcg    = 4;                   % ECG disturbance threshold (applies to both ECG scores)
thrAcc    = 4;
padAccSec = 2;
padEcgSec = 2;
% Which regions to filter
targetMode  = 'orange' ;       % 'orange' = all aligned ACC motion, 'red' = motion AND ECG disturbance
includeEcgOnly = true;           % 'orange' mode only: also filter ECG disturbances with no ACC motion
% Filters
medWinSec    = [0.2 0.6];        % cascaded median filters used to estimate baseline wander (s)
mildLowHz    = 40;               % mild tier low-pass corner (Hz)
heavyBandHz  = [5 25];           % heavy tier band-pass corners (Hz)
filtOrder    = 3;                % Butterworth order (band-pass doubles it)
taperSec     = 0.25;             % cross-fade length at region edges (s)
% Optional wavelet shrinkage inside the heavy tier (OFF by default).
% In tests on S002 it did NOT beat the 5-25 Hz band-pass alone; try it and compare.
useWavelet   = false;
waveName     = 'sym4';
waveLevel    = 5;                % stationary wavelet transform levels (100 Hz: d1 25-50 Hz ... d4 3-6 Hz)
waveK        = 1.0;              % soft-threshold = waveK x local noise level (lower = gentler)
waveMarginSec = 2;               % context added each side of every heavy region
% Severe-artifact handling (filtering can't recover these)
blankFactor = 6;                 % severe if window's ECG deviation > blankFactor x typical
blankPadSec = 1;                 % widen severe windows by this much each side
blankSevere = false;             % true = set severe windows to NaN in the output
% Rhythm-constrained beat selection
useRhythm      = true;
rhythm.minScore = 0.30;          % candidate must score above this to be worth including
rhythm.alpha    = 0.5;           % penalty weight for RR deviating from the local rhythm
rhythm.beta     = 1.5;           % penalty per missed beat
rhythm.kmax     = 4;             % allow up to kmax-1 consecutive missed beats
rhythm.sigFrac  = 0.10;          % expected RR spread as a fraction of local RR
rhythm.maxDev   = 3.5;           % reject RR further than this many spreads from the rhythm
% Display / output
viewStart = [];                  % [] = centre the view on the strongest flagged region
viewLen   = 30;
saveFiltered = false;            % true = write <name>_filtered.csv (original file is never modified)
%% ---- Load original data -------------------------------------------------
if ~isfile(filename)
   [f, p] = uigetfile({'*.csv', 'CSV files'}, 'Select 100 Hz DREAMT file');
   if isequal(f, 0), error('No file selected.'); end
   filename = fullfile(p, f);
end
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
opts.SelectedVariableNames = {'TIMESTAMP', 'ECG', 'ACC_X', 'ACC_Y', 'ACC_Z'};
T = readtable(filename, opts);
t    = T.TIMESTAMP(:);
ecg  = T.ECG(:);
accX = T.ACC_X(:);  accY = T.ACC_Y(:);  accZ = T.ACC_Z(:);
good = ~isnan(t);
[t, order] = sort(t(good));
ecg  = ecg(good);  ecg  = ecg(order);
accX = accX(good); accX = accX(order);
accY = accY(good); accY = accY(order);
accZ = accZ(good); accZ = accZ(order);
fs = 1 / median(diff(t));
N  = numel(t);
fprintf('Rows: %d | %.2f hours | fs = %.1f Hz\n', N, (t(end) - t(1)) / 3600, fs);
ecg  = fillmissing(ecg,  'linear');
accX = fillmissing(accX, 'previous');
accY = fillmissing(accY, 'previous');
accZ = fillmissing(accZ, 'previous');
if max(abs(ecg)) < 0.01, ecg = ecg * 1e6; end             % V -> uV
accX = accX / 64;  accY = accY / 64;  accZ = accZ / 64;   % 1/64 g -> g
tp = t - t(1);
%% ---- 1-2) Window scores and ACC->ECG delay ------------------------------
w = round(winSec * fs);
n = floor(N / w);
twp = (0:n-1)' * winSec;
E  = reshape(ecg(1:n*w),  w, n);
Ax = reshape(accX(1:n*w), w, n);
Ay = reshape(accY(1:n*w), w, n);
Az = reshape(accZ(1:n*w), w, n);
ecgMetric = max(abs(E - median(E, 1)), [], 1)';                           % amplitude disturbance
accMetric = sum(sqrt(diff(Ax, 1, 1).^2 + diff(Ay, 1, 1).^2 + diff(Az, 1, 1).^2), 1)';
% High-frequency (25-45 Hz) noise: catches EMG-type bursts whose amplitude is
% similar to the R peaks and so would be missed by the amplitude score alone.
[zh_, ph_, kh_] = butter(4, [25 45] / (fs / 2), 'bandpass');
[sosHF, gHF] = zp2sos(zh_, ph_, kh_);
ecgHF = filtfilt(sosHF, gHF, ecg);
H = reshape(ecgHF(1:n*w), w, n);
hfMetric = sqrt(mean(H.^2, 1))';
ze = robustZ(ecgMetric);
zh = robustZ(hfMetric);
zEcgAny = max(ze, zh);                      % disturbed by either measure
za = robustZ(accMetric);
maxLagW = round(maxLagSec / winSec);
[lags, cc] = lagCorr(ze, za, maxLagW);
[rBest, ib] = max(cc);
lagEst = lags(ib) * winSec;
if isempty(lagOverrideSec), lagSec = lagEst; else, lagSec = lagOverrideSec; end
lagW = round(lagSec / winSec);
fprintf ('Delay (positive = ACC after ECG): estimated %+.1f s (r = %.2f), using %+.1f s\n', lagEst, rBest, lagSec);
%% ---- 3) Flag regions ----------------------------------------------------
zaAl = nan(n, 1);
if lagW >= 0, zaAl(1:n-lagW) = za(1+lagW:n); else, zaAl(1-lagW:n) = za(1:n+lagW); end
kAcc = max(0, round(padAccSec / winSec));
kEcg = max(0, round(padEcgSec / winSec));
motionDil = movmax(double(zaAl > thrAcc), 2*kAcc + 1) > 0;
ecgDil    = movmax(double(zEcgAny > thrEcg), 2*kEcg + 1) > 0;
switch lower(targetMode)
   case 'orange'
       if includeEcgOnly, flagged = motionDil | ecgDil; else, flagged = motionDil; end
   case 'red'
       flagged = motionDil & ecgDil;
   otherwise
       error('targetMode must be ''orange'' or ''red''.');
end
heavyMask = flagged & ecgDil;               % heavy tier only where the ECG is actually disturbed
mildMask  = flagged;
fprintf('Flagged (mild tier):  %d windows (%.2f%% of recording)\n', sum(mildMask),  100 * mean(mildMask));
fprintf('Heavy tier:           %d windows (%.2f%% of recording)\n', sum(heavyMask), 100 * mean(heavyMask));
fprintf('Aligned ACC motion where ECG looks clean: %d windows (mild tier only)\n', sum(motionDil & ~ecgDil));
%% ---- 4) Tiered filtering inside flagged regions -------------------------
k1 = 2 * round(medWinSec(1) * fs / 2) + 1;      % odd kernel lengths
k2 = 2 * round(medWinSec(2) * fs / 2) + 1;
baseline = medfilt1(medfilt1(ecg, k1, 'truncate'), k2, 'truncate');
ecgNB = ecg - baseline;                          % baseline wander removed
[zl, pl, kl] = butter(filtOrder, mildLowHz / (fs / 2), 'low');
[sosL, gL]   = zp2sos(zl, pl, kl);
mildSig = filtfilt(sosL, gL, ecgNB);             % zero-phase: no R-peak shift
[zb, pb, kb] = butter(filtOrder, heavyBandHz / (fs / 2), 'bandpass');
[sosB, gB]   = zp2sos(zb, pb, kb);
heavySig = filtfilt(sosB, gB, ecgNB);
if useWavelet
   if isempty(which('swt')) || ~license('test', 'Wavelet_Toolbox')
       warning('Wavelet Toolbox not available - skipping wavelet shrinkage.');
   else
       fprintf('Applying wavelet shrinkage inside %d heavy-tier regions...\n', ...
               size(maskToIntervals(heavyMask, winSec), 1));
       heavySig = waveletShrinkRegions(heavySig, heavyMask, winSec, fs, waveName, waveLevel, waveK, waveMarginSec);
   end
end
taperN = max(1, round(taperSec * fs));
w1 = movmean(double(winToSamp(mildMask,  w, N)), taperN);     % 1 inside flagged regions
w2 = min(movmean(double(winToSamp(heavyMask, w, N)), taperN), w1);
ecgOut = ecg + w1 .* (mildSig - ecg) + w2 .* (heavySig - mildSig);
%% ---- 5) Severe windows --------------------------------------------------
kBlank = max(0, round(blankPadSec / winSec));
severeWin = flagged & (ecgMetric > blankFactor * median(ecgMetric));
severeDil = movmax(double(severeWin), 2*kBlank + 1) > 0;
fprintf('Severe-artifact windows: %d (+padding: %d)  | blanking %s\n', ...
       sum(severeWin), sum(severeDil), ternary(blankSevere, 'ON', 'OFF'));
if blankSevere && any(severeDil)
   ecgOut(winToSamp(severeDil, w, N)) = NaN;
end
%% ---- 6) Beat detection: naive vs rhythm-constrained ---------------------
% Naive detector has NO internal band-pass (so linear pre-filtering can't be
% "undone" by the detector): prominence threshold from clean-window R amplitude.
promThr = 0.5 * median(ecgMetric(~ecgDil & ~motionDil));
pkRaw  = detectR(ecg,    fs, promThr);
pkFilt = detectR(ecgOut, fs, promThr);
sRaw   = hrStats(pkRaw,  t);
sFilt  = hrStats(pkFilt, t);
R = struct('idx', pkFilt(:), 'nStretches', 0, 'nUnresolved', 0);
if useRhythm
   fprintf('\nRhythm selection: %d uncertain stretches, %d unresolved (kept naive detections)\n', ...
       R.nStretches, R.nUnresolved);
end
sRhy = hrStats(R.idx, t);
% QRS-to-noise inside flagged regions, measured at the same beat positions in all signals
flaggedSamp = winToSamp(mildMask, w, N);
refPk = R.idx;
qRaw  = qrsToNoise(ecg,    refPk, flaggedSamp);
qFilt = qrsToNoise(ecgOut, refPk, flaggedSamp);
fprintf('\n%-38s %10s %10s %12s\n', '', 'Original', 'Filtered', '+ Rhythm');
fprintf('%-38s %10d %10d %12d\n',     'R peaks',                              sRaw.beats, sFilt.beats, sRhy.beats);
fprintf('%-38s %10.1f %10.1f %12.1f\n','Average HR (bpm)',                    sRaw.avg,   sFilt.avg,   sRhy.avg);
fprintf('%-38s %10.1f %10.1f %12.1f\n','Min instantaneous HR (bpm)',          sRaw.min,   sFilt.min,   sRhy.min);
fprintf('%-38s %10.1f %10.1f %12.1f\n','Max instantaneous HR (bpm)',          sRaw.max,   sFilt.max,   sRhy.max);
fprintf('%-38s %10d %10d %12d\n',     'Implausible beats (HR <40 or >140)',   sRaw.bad,   sFilt.bad,   sRhy.bad);
fprintf('%-38s %10.1f %10.1f\n',       'QRS-to-noise in flagged regions',     qRaw,       qFilt);
if useRhythm
   fprintf('\nRhythm selection: %d uncertain stretches, %d unresolved (kept naive detections), %d interpolated beats\n', ...
           R.nStretches, R.nUnresolved);
   fprintf('Interpolated beats are timing estimates from the rhythm, NOT detected QRS complexes -\n');
   fprintf('exclude them (result.beatIsInterp) from beat-to-beat variability analysis.\n');
end
%% ---- Optional: save (original CSV is never touched) ---------------------
if saveFiltered
   [pth, nm] = fileparts(filename);
   outCsv = fullfile(pth, [nm '_filtered.csv']);
   writetable(table(t, ecg, ecgOut, 'VariableNames', {'TIMESTAMP', 'ECG_original_uV', 'ECG_filtered_uV'}), outCsv);
   fprintf('Wrote %s\n', outCsv);
end
%% ---- Plot: original vs filtered + aligned accelerometer + disturbance scores ----
ivMild   = maskToIntervals(mildMask,   winSec);
ivHeavy  = maskToIntervals(heavyMask,  winSec);
ivSevere = maskToIntervals(severeDil,  winSec);
if isempty(viewStart)
   if any(heavyMask)
       zt = zEcgAny; zt(~heavyMask) = -inf;
       [~, jMax] = max(zt);
       viewStart = max(0, twp(jMax) - viewLen / 2);
   else
       viewStart = 0;
   end
end
xl = [0, 0 + viewLen];
tAcc = tp - lagSec;                          % ACC shifted onto the ECG time axis
figure('Name', 'Original vs filtered ECG', 'Color', 'w', 'Position', [100 30 1250 1100]);
% -- Panel 1: Original ECG
ax1 = subplot(4, 1, 1);
hE1 = plot(tp, ecg, 'b', 'LineWidth', 0.8); hold on;
hP1 = plot(tp(pkRaw), ecg(pkRaw), 'ko', 'MarkerSize', 7, 'LineWidth', 1.2);
hM1 = shade(ax1, ivMild,   [1.00 0.80 0.35], 0.35);
hH1 = shade(ax1, ivHeavy,  [0.90 0.20 0.20], 0.25);
hS1 = shade(ax1, ivSevere, [0.35 0.35 0.35], 0.30);
hold off; grid on;
ylabel('ECG (\muV)');
title(sprintf('Original ECG  |  %d R peaks, %d implausible', sRaw.beats, sRaw.bad));
[h1, n1] = legendEntries({hE1, hP1, hM1, hH1, hS1}, ...
   {'Original ECG', 'Detected R peaks', 'Mild tier', 'Heavy tier', 'Severe windows'});
legend(h1, n1, 'Location', 'northeast');
% -- Panel 2: Filtered ECG
ax2 = subplot(4, 1, 2);
hE2 = plot(tp, ecgOut, 'Color', [0 0.5 0.1], 'LineWidth', 0.8); hold on;
hP2 = plot(tp(R.idx), ecgOut(R.idx), 'ko', 'MarkerSize', 7, 'LineWidth', 1.2);   % black circles
hM2 = shade(ax2, ivMild,   [1.00 0.80 0.35], 0.35);
hH2 = shade(ax2, ivHeavy,  [0.90 0.20 0.20], 0.25);
hS2 = shade(ax2, ivSevere, [0.35 0.35 0.35], 0.30);
hold off; grid on;
ylabel('ECG (\muV)');
title(sprintf('Filtered%s%s  |  %d R peaks, %d implausible', ...
   ternary(useWavelet, ' + wavelet', ''), ternary(useRhythm, ' + rhythm selection', ''), sRhy.beats, sRhy.bad));
[h2, n2] = legendEntries({hE2, hP2, hM2, hH2, hS2}, ...
   {'Filtered ECG', 'Detected R peaks', 'Mild tier', 'Heavy tier', 'Severe windows'});
legend(h2, n2, 'Location', 'northeast');
% -- Panel 3: Accelerometer (shifted to line up with ECG)
ax3 = subplot(4, 1, 3);
hX3 = plot(tAcc, accX, 'r', 'LineWidth', 1); hold on;
hY3 = plot(tAcc, accY, 'g', 'LineWidth', 1);
hZ3 = plot(tAcc, accZ, 'b', 'LineWidth', 1);
shade(ax3, ivMild,   [1.00 0.80 0.35], 0.35);
shade(ax3, ivHeavy,  [0.90 0.20 0.20], 0.25);
shade(ax3, ivSevere, [0.35 0.35 0.35], 0.30);
hold off; grid on;
ylabel('Acceleration (g)');
title(sprintf('Accelerometer (shifted by %+.1f s to line up with ECG)', lagSec));
legend([hX3 hY3 hZ3], {'ACC\_X', 'ACC\_Y', 'ACC\_Z'}, 'Location', 'northeast');
% -- Panel 4: Disturbance scores per window
ax4 = subplot(4, 1, 4);
hZe = stairs(twp, ze, 'b', 'LineWidth', 1); hold on;
hZh = stairs(twp, zh, 'Color', [0.30 0.75 0.93], 'LineWidth', 0.8);
hZa = stairs(twp, zaAl, 'Color', [0.85 0.45 0.0], 'LineWidth', 1);
plot([0 tp(end)], [thrEcg thrEcg], 'b--');
plot([0 tp(end)], [thrAcc thrAcc], '--', 'Color', [0.85 0.45 0.0]);
shade(ax4, ivMild,   [1.00 0.80 0.35], 0.35);
shade(ax4, ivHeavy,  [0.90 0.20 0.20], 0.25);
shade(ax4, ivSevere, [0.35 0.35 0.35], 0.30);
hold off; grid on;
ylabel('Robust z-score');
xlabel('Time since start of recording (s)');
title('Disturbance scores per window (dashed = thresholds)');
legend([hZe hZh hZa], {'ECG amplitude disturbance', 'ECG high-frequency noise', 'ACC motion (aligned)'}, ...
   'Location', 'northeast');
ylim(ax4, [-2, max([ze; zh; zaAl; thrEcg + 2; thrAcc + 2], [], 'omitnan') * 1.05]);
% -- Linking: x shared by all four panels, y shared only by the two ECG panels
linkaxes([ax1 ax2 ax3 ax4], 'x');
yLink = linkprop([ax1 ax2], 'YLim');
setappdata(gcf, 'ecgYLink', yLink);          % keep the link object alive
xlim(ax1, xl);
ax1.XAxis.Exponent = 0;  ax2.XAxis.Exponent = 0;
ax3.XAxis.Exponent = 0;  ax4.XAxis.Exponent = 0;
% Fit the y-ranges to the initial view
i1 = tp >= xl(1) & tp <= xl(2);
setYFit(ax1, ecg(i1));                       % both ECG panels share the original's y-range
i3 = tAcc >= xl(1) & tAcc <= xl(2);
setYFit(ax3, [accX(i3); accY(i3); accZ(i3)]);
%% ---- Results struct for further work ------------------------------------
result = struct('t', t, 'ecgOriginal', ecg, 'ecgFiltered', ecgOut, 'fs', fs, 'lagSec', lagSec, ...
   'beatIdx', R.idx, ...
   'mildRegionsSec', ivMild, 'heavyRegionsSec', ivHeavy, 'severeRegionsSec', ivSevere);
