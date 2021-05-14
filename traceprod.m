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
function y = traceprod(A, B, option)
%
% traceprod(A,B) = trace(A*B).  The trace can be computed more efficiently depending on the value of option.
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
  if strcmp(option, 'Isotropic')
    y = A*trace(B);
  elseif strcmp(option, 'Diag')
    y = A'*diag(B);
  elseif strcmp(option, 'Full')
    y = trace(A*B);
  end
end
