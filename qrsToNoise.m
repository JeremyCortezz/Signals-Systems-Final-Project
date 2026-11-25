function q = qrsToNoise(x, refPk, mask)
   % Median R amplitude at the reference beats inside the mask, divided by the
   % robust noise level (1.4826 x MAD) of the signal inside the mask.
   x = fillmissing(x, 'constant', 0);
   pk = refPk(mask(refPk));
   if numel(pk) < 3 || nnz(mask) < 100, q = NaN; return; end
   amp = zeros(numel(pk), 1);
   for k = 1:numel(pk)
       i1 = max(1, pk(k) - 4);  i2 = min(numel(x), pk(k) + 4);
       amp(k) = max(x(i1:i2));
   end
   xm = x(mask);
   q = median(amp) / (1.4826 * median(abs(xm - median(xm))) + eps);
end
