function y = waveletShrinkRegions(y, mask, winSec, fs, wname, L, k, marginSec)
   % Stationary-wavelet soft shrinkage with a TIME-LOCAL noise estimate
   % (moving median of |detail| / 0.6745), applied region by region so the
   % threshold adapts to each noise burst. Drops the level-L detail and the
   % approximation (< ~1.6 Hz), keeping levels 1..L-1.
   y0 = y;
   N  = numel(y);
   iv = maskToIntervals(mask, winSec);
   pad = 2^L;
   m   = round(marginSec * fs);
   for r = 1:size(iv, 1)
       i1 = max(1, round(iv(r, 1) * fs) + 1 - m);
       i2 = min(N, round(iv(r, 2) * fs) + m);
       seg = y0(i1:i2);
       Ls  = numel(seg);
       Lp  = ceil(Ls / pad) * pad;
       if Lp > Ls, segp = [seg; seg(end:-1:end - (Lp - Ls) + 1)]; else, segp = seg; end
       [swa, swd] = swt(segp, L, wname);
       for l = 1:L
           if l == L
               swd(l, :) = 0;
           else
               d = swd(l, :);
               sg = movmedian(abs(d), 101) / 0.6745;
               swd(l, :) = sign(d) .* max(abs(d) - k * sg, 0);
           end
       end
       swa(L, :) = 0;
       rec = iswt(swa, swd, wname);
       y(i1:i2) = rec(1:Ls);
   end
end
