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
function model = initialize_parameters(I, J)
%
% Generate parameters of MLDS model
%
% ---Inputs---
% I:  observation dimensionality, a vector of positive integers
% J:  latent dimensionality, a vector of positive integers
%
% ---Outputs---
% model.u0:  vec(u0), where u0 is the expectation tensor of Z_1
% model.Q0:  mat(Q0), where Q0 is the covariance tensor of Z_1
% model.Q:  mat(Q), where Q is the covariance tensor of Z_n+1 | Z_n
% model.R:  mat(R), where R is the covariance tensor of X_n | Z_n
% model.A:  cell array of projection matrices of the transition multilinear operator A
% model.C:  cell array of projection matrices of the projection multilinear operator C
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
  M = numel(J);
  I = reshape(I,1,M);
  J = reshape(J,1,M);
  prodI = prod(I);
  prodJ = prod(J);
  ItimesJ = I .* J;
  
  model.mu0 = vec2ten(randn(prodJ,1), J);
  	
  model.Q0 = mat2ten(eye(prodJ), [J J]);
  model.Q = mat2ten(eye(prodJ), [J J]);
  model.R = mat2ten(eye(prodI), [I I]);
  model.A = initialize_multilinear_operator(J,J);
  model.C = initialize_multilinear_operator(I,J);
end

%-----------------------------------------------------
function C = initialize_multilinear_operator(I, J)
  prodI = prod(I);
  prodJ = prod(J);
  M = numel(I);
  C = cell(M,1);
  for m = 1:M
    rC = randn(I(m));
    while rank(rC) < I(m)
      rC = randn(I(m));
    end
    [U S V] = svd(rC);
    C{m} = U(:,1:J(m));
  end
  if M == 1
    C = C{1};
  end
end
