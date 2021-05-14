% This source code is (c) Copyright by Lei Li.
% All rights preserved.
%
% Permission is granted to use it for non-profit purposes,
% including research and teaching. For-profit use requires
% the express consent of the author (leili@cs.cmu.edu).

function [Y] = linear_interp(X, W)
% use linear interpolation to fill in the missing value
% W is the indication matrix, 1=observed, 0=missing.
% use the matlab internal interpolation method
% modified by leili (2009-5-5)
% modified by leili (2010-1-21)
% modified by leili (2010-4-5),
%   modified the extrapolation to be nearest neighbor
% modified by leili (2013-1-4),
%   if the extrapolation part is empty, then no need to extrapolation

M = size(X, 1);
N = size(X, 2);
idx = 1:N;
Y = X;
for i = 1:M
  obs = W(i, :) ~= 0;
  if sum(obs) > 1
    yy = interp1(idx(obs), X(i, obs), idx(~obs), 'linear');
    Y(i, ~obs) = yy;
    % using neareast neighbor for extrapolation
    idyy = idx(isnan(Y(i, :)));
    if ~isempty(idyy)
      yy = interp1(idx(obs), X(i, obs), idyy, 'nearest', 'extrap');
      Y(i, idyy) = yy;
    end
  else
    Y(i, ~obs) = 0;
  end
end


% $$$ M = size(X,2);
% $$$ N = size(X,1);
% $$$
% $$$ Y = X;
% $$$ for i = 1:M
% $$$     s = findFirst(W(:, i), 1);
% $$$     t = findFirst(W(:, i), s+1);
% $$$     if (t > 2)
% $$$         k = (X(t, i) - X(s,i)) / (t-s);
% $$$         for j = 1:(t-1)
% $$$             Y(j, i) = X(s, i) + (j - s) * k;
% $$$         end
% $$$     end
% $$$     while (t < N)
% $$$         olds = s;
% $$$         s = t;
% $$$         t = findFirst(W(:,i), t+1);
% $$$         if (t <= N)
% $$$             k = (X(t, i) - X(s,i)) / (t-s);
% $$$             for j = (s+1):(t-1)
% $$$                 Y(j, i) = X(s, i) + (j - s) * k;
% $$$             end
% $$$         else
% $$$             t = s;
% $$$             s = olds;
% $$$             k = (X(t, i) - X(s,i)) / (t-s);
% $$$             for j = (t+1):N
% $$$                 Y(j, i) = X(s, i) + (j - s) * k;
% $$$             end
% $$$             t = N+1;
% $$$         end
% $$$     end
% $$$ end
% $$$
% $$$ function [idx] = findFirst(f, a)
% $$$ while ((a <= size(f, 1)) && (f(a) == 0))
% $$$     a = a+1;
% $$$ end
% $$$ idx = a;
% $$$
% $$$
