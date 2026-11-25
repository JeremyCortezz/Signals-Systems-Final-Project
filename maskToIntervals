function iv = maskToIntervals(mask, winSec)
   d = diff([0; double(mask(:)); 0]);
   s = find(d == 1);
   e = find(d == -1) - 1;
   iv = [(s - 1) * winSec, e * winSec];
end
