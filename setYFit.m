function setYFit(ax, y)
   if isempty(y), return; end
   lo = min(y);  hi = max(y);
   pad = 0.1 * (hi - lo);
   if pad == 0, pad = 1e-3; end
   ylim(ax, [lo - pad, hi + pad]);
end
