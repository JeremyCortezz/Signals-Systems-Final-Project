function pk = detectR(x, fs, promThr)
   % Prominence-based detector on the signal as given (no internal band-pass).
   x = fillmissing(x, 'constant', 0);
   [~, pk] = findpeaks(x, 'MinPeakProminence', promThr, ...
                          'MinPeakDistance', round(0.40 * fs));
end
