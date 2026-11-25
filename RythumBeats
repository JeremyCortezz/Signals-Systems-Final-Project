function R = rhythmBeats(y, fs, w, naivePk, uncertainWin, par)
% Inside each "uncertain" stretch, pick the subset of QRS-shaped candidates that
% best fits the heart rhythm set by trusted beats on either side (dynamic
% programming). Outside those stretches the naive detections are kept.
N = numel(y);
naivePk = naivePk(:);
uncSamp = winToSamp(uncertainWin, w, N);
trusted = naivePk(~uncSamp(naivePk));
R = struct('idx', naivePk, 'nStretches', 0, 'nUnresolved', 0);
if numel(trusted) < 20 || ~any(uncertainWin), return; end
% QRS template and typical amplitude from trusted beats
sel = trusted(trusted > 12 & trusted < N - 12);
if numel(sel) > 2000, sel = sel(round(linspace(1, numel(sel), 2000))); end
S = zeros(21, numel(sel));
A = zeros(numel(sel), 1);
for q = 1:numel(sel)
   S(:, q) = y(sel(q) - 10 : sel(q) + 10);
   A(q) = max(y(sel(q) - 4 : sel(q) + 4));
end
tmpl = median(S, 2);
Aref = median(A);
tz = tmpl - mean(tmpl);  tz = tz / norm(tz);
m  = numel(tz);
% Shape score: normalised cross-correlation with the template
num = conv(y, flipud(tz), 'same');
s1  = conv(y,    ones(m, 1), 'same');
s2  = conv(y.^2, ones(m, 1), 'same');
ncc = num ./ sqrt(max(s2 - s1.^2 / m, 1e-9));
[~, cand] = findpeaks(ncc, 'MinPeakHeight', 0.3, 'MinPeakDistance', round(0.25 * fs));
iv = maskToIntervals(uncertainWin, w / fs);
R.nStretches = size(iv, 1);
keepIdx = trusted;
for r = 1:size(iv, 1)
   i1 = max(1, round(iv(r, 1) * fs) + 1);
   i2 = min(N, round(iv(r, 2) * fs));
   a0 = find(trusted < i1, 1, 'last');
   a1 = find(trusted > i2, 1, 'first');
   inStretch = naivePk(naivePk >= i1 & naivePk <= i2);
   resolved = false;
   if ~isempty(a0) && ~isempty(a1)
       tb = trusted(max(1, a0 - 8) : a0);
       ta = trusted(a1 : min(numel(trusted), a1 + 8));
       rrs = [diff(tb); diff(ta)] / fs;
       rrs = rrs(rrs > 0.4 & rrs < 1.6);
       c = cand(cand > trusted(a0) & cand < trusted(a1));
       if numel(rrs) >= 3 && numel(c) >= 1 && numel(c) <= 1500
           rrref = median(rrs);
           amp = zeros(numel(c), 1);
           for q = 1:numel(c)
               amp(q) = max(y(max(1, c(q) - 3) : min(N, c(q) + 3)));
           end
           score = ncc(c) .* min(1, amp / Aref);
           [ok, selB] = selectBeatsDP(c, score, trusted(a0), trusted(a1), rrref, fs, par);
           if ok
               resolved = true;
               keepIdx  = [keepIdx; selB(:)];          %#ok<AGROW>
           end
       end
   end
   if ~resolved
       R.nUnresolved = R.nUnresolved + 1;
       keepIdx = [keepIdx; inStretch];                 %#ok<AGROW>
   end
end
R.idx = sort(keepIdx);
R.idx = min(max(R.idx, 1), N);
end
