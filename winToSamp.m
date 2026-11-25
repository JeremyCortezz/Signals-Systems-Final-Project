function s = winToSamp(mask, w, N)
   s = repelem(logical(mask(:)), w);
   s = [s; false(N - numel(s), 1)];
end
