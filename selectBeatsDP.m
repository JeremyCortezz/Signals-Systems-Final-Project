function [ok, sel] = selectBeatsDP(cand, score, a0, a1, rrref, fs, par)
% Dynamic programming over candidate beats between two trusted anchor beats.
% gain  = shape score - threshold (reward for including a QRS-like candidate)
% cost  = alpha*((RR/k - RRref)/sigma)^2 + beta*(k-1)   (k-1 = skipped beats between)
sig   = par.sigFrac * rrref;
nodes = [a0; cand(:); a1];
gain  = [0; score(:) - par.minScore; 0];
M = numel(nodes);
best = -inf(M, 1);  prev = zeros(M, 1);
best(1) = 0;
for j = 2:M
   i  = (1:j-1)';
   rr = (nodes(j) - nodes(i)) / fs;
   cost = inf(j-1, 1);
   for k = 1:par.kmax
       dev = rr / k - rrref;
       c = par.alpha * (dev / sig).^2 + par.beta * (k - 1);
       c(abs(dev) > par.maxDev * sig) = inf;
       cost = min(cost, c);
   end
   cost(rr < 0.30 | ~isfinite(best(i))) = inf;
   tot = best(i) + gain(j) - cost;
   [v, ib] = max(tot);
   if isfinite(v)
       best(j) = v;  prev(j) = ib;
   end
end
ok = isfinite(best(M));
sel = [];
if ~ok, return; end
j = prev(M);                     % walk back from the end anchor, skip both anchors
while j > 1
   sel(end+1, 1) = nodes(j);    %#ok<AGROW>
   j = prev(j);
end
sel = flipud(sel);
end
