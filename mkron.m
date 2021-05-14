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
function mkronA = mkron(A)
%
% if A is a cell array, return A{M} * A{M-1} * ... * A{1}, where * is the Kronecker matrix product.  If A is already a matrix, then return A.
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
if ~strcmp(class(A),'cell')
  mkronA = A;
elseif numel(A) == 1
  mkronA = A{1};
else
  M = numel(A);
  mkronA = kron(A{M}, mkron(subcell(A,1:M-1)));
end
