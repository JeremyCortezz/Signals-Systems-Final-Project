function s = hrStats(pk, t)
   s = struct('beats', numel(pk), 'avg', NaN, 'min', NaN, 'max', NaN, 'bad', 0);
   if numel(pk) < 2, return; end
   rr  = diff(t(pk));
   bpm = 60 ./ rr;
   s.avg = 60 / mean(rr);
   s.min = min(bpm);
   s.max = max(bpm);
   s.bad = sum(bpm < 40 | bpm > 140);
end
