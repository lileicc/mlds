% This source code is (c) Copyright by Lei Li, Mark Rogers.
% All rights preserved.
%
% Permission is granted to use it for non-profit purposes,
% including research and teaching. For-profit use requires
% the express consent of the author (leili@cs.berkeley.edu).
%
% Details in the following paper:
%   Mark Rogers, Lei Li and Stuart J. Russell (2013),
%     "Multilinear Dynamical Systems for Tensor Time Series",    
%     In Advances in Neural Information Processing Systems 26. 
%
function v = ten2vec(T)
%
% ten2vec is an overloaded function such that
% if T is a tensor:
%   - return vec(T)
% if T is a cell array of tensors, e.g., a tensor time series:
%   - return a matrix v such that v(:,n) = vec(T{n})
% if T is a struct, e.g., the MLDS model parameters:
%   - return the same struct but with each parameter vectorized or matricized. To preserve the projection matrices A{m} and C{m}, attributes v.cellA and v.cellC are created.
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
switch class(T)
case 'double'
  v = vec(T);
case 'cell'
  p = prod(size(T{1}));
  N = numel(T);
  v = zeros(p, N);
  for n = 1:N
    v(:, n) = ten2vec(T{n});
  end
case 'struct'
  v.mu0 = ten2vec(T.mu0);
  v.Q0 = ten2mat(T.Q0);
  v.Q = ten2mat(T.Q);
  v.R = ten2mat(T.R);
  v.cellA = T.A;
  v.cellC = T.C;
  v.A = mkron(T.A);
  v.C = mkron(T.C);
end
