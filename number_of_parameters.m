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
function p = number_of_parameters(I, J, Type)
%
% count the number of MLDS parameters
%
% ---Inputs---
% I:  observation dimensionality, vector of positive integers
% J:  latent dimensionality, vector of positive integers
% Type.Q0:  either 'Isotropic', 'Diag', or 'Full'
% Type.Q:  either 'Isotropic', 'Diag', or 'Full'
% Type.R:  either 'Isotropic', 'Diag', or 'Full'
%
% ---Outputs---
% p:  number of parameters
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
  M = numel(I);
  I = reshape(I,1,M);
  J = reshape(J,1,M);
  prodI = prod(I);
  prodJ = prod(J);
  p = count_covariance_parameters(J, Type.Q0) ...
    + count_covariance_parameters(J, Type.Q) ...
    + count_covariance_parameters(I, Type.R) ...
    + sum(J .* J) + sum(I .* J);
end

%---------------------------------------------------

function p = count_covariance_parameters(I, Type)
  p = 0;
  switch Type
  case 'Isotropic'
    p = p + 1;
  case 'Diag'
    p = p + prod(I);
  case 'Full'
    p = p + prod(I)^2;
  end
end
