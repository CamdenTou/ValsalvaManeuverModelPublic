function [y, train] = peaks2(x, nbd, thresh)
% PEAKS(x, nbd)
% Find all the peaks in vector x that are local maxes and at least
% as high as any other point within nbd elements.  Returns a vector
% of indices into x.  If a peak is broad, uses the lowest index of
% the several possible.
%
% PEAKS(x, nbd, thresh)
% As above, but a peak must additionally be at least as big as thresh.
%
% [y,train] = PEAKS( ... )
% If there are two return values, the second is a "spike train" the same
% length as x that has values from x at x's peaks and zeroes elsewhere.

x = x(:);
n = length(x);
% endpoints can be peaks too
if (x(1) >= x), p0 = x(1) > x(n);
else            p0 = (x(1) > x(2)); end
p1 = (x(n-1) < x(n));

% p is the initial guess at the peak locations
p = find([p0; (x(1:n-2) < x(2:n-1)) & (x(2:n-1) >= x(3:n)); p1]);
if (nargin > 2), p = p(find(x(p) >= thresh)); end

% Now test points in p for being greater than neighbors.
v = x(p);
y = [];
for i = 1:length(p)
   if (v(i) >= v(find( abs(p - p(i)) <= nbd ))),
      if (v(i) >= x(max(1,p(i)-nbd) : min(n,p(i)+nbd))),
         y = [y; p(i)];
      end
   end
end

if (nargout > 1),
   train = zeros(n, 1);
   train(y) = x(y);
end