function [lags, cc] = lagCorr(a, b, maxLag)
   lags = -maxLag:maxLag;
   cc = nan(size(lags));
   n = numel(a);
   for k = 1:numel(lags)
       L = lags(k);
       if L >= 0, x = a(1:n-L);   y = b(1+L:n);
       else,      x = a(1-L:n);   y = b(1:n+L); end
       if numel(x) > 10
           c = corrcoef(x, y);
           cc(k) = c(1, 2);
       end
   end
end
