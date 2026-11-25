function z = robustZ(x)
   x = log1p(x - min(x));
   m = median(x);
   s = 1.4826 * median(abs(x - m)) + eps;
   z = (x - m) / s;
end
