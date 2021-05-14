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
function T = vec2ten(v, varargin)
%
% vec2ten is an overloaded function such that
% if v is a vector:
%   - return T with size I = varargin{1} such that vec(T) = v,
% if v is a matrix, e.g., a vectorized tensor time series:
%   - return a cell array T such that vec(T{n}) = v(:,n), where size(T{n}) = I = varargin{1}
% if v is a struct, e.g., the vectorized MLDS model parameters:
%   - return the same struct but with each parameters unvectorized or unmatricized.
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
switch class(v)
case {'double' 'logical'}
  I = varargin{1};
  if numel(v) == prod(I)
    T = zeros([I 1]);
    T(:) = v;
  else
    N = size(v,2);
    T = cell(N,1);
    for n = 1:N
      T{n} = vec2ten(v(:,n), I);
    end
  end
case 'struct'
  switch class(v.cellA)
  case 'double'
    I = size(v.cellC,1);
    J = size(v.cellC,2);
  case 'cell'
    M = numel(v.cellA);
    I = zeros(1,M);
    J = zeros(1,M);
    for m = 1:M
      I(m) = size(v.cellC{m},1);
      J(m) = size(v.cellC{m},2);
    end
  end
  T.mu0 = vec2ten(v.mu0, J);
  T.Q0 = mat2ten(v.Q0, [J J]);
  T.Q = mat2ten(v.Q, [J J]);
  T.R = mat2ten(v.R, [I I]);
  T.A = v.cellA;
  T.C = v.cellC;
end
