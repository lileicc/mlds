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
function A = ten2mat(T)
%
% A = matricize(T)
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
I = size(T);
M = numel(I);
if mod(M,2) == 1
  M = M + 1;
  I(M) = 1;
end
A = zeros(prod(I(1:(M/2))), prod(I([1:(M/2)]+(M/2))));
A(:) = T(:);
