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
function [model diagnostics] = learn_mlds(X, varargin)
%
% fit MLDS to data X
%
% ---Inputs---
% X:  tensor time series represented by an [N x 1] cell array. size(X{n}) = I for all n, where I is the observation dimensionality, a vector of positive integers.
%
% ---Optional inputs---
% 'J', followed a vector of positive integers J indicating the latent dimensionality.  J = I by default.
% 'W', followed by a cell array W of tensors of length N such that each entry of W{n} is 1 if missing and 0 otherwise.  W{n} = zeros(size(I)) for all n by default.
% 'MaxIter', followed by a positive integer MaxIter indicating the maximum number of EM iterations.  MaxIter = 10 by default.
% 'Epsilon', followed by a scalar Epsilon indicating the convergence bound for the log-likelihood. Epsilon = 1e-5 by default.
% 'Type', followed by a a struct Type with attributes indicating the covariance tensors such that
%  - Type.Q0 is equal to one of 'Isotropic', 'Diag', 'Full'.  Type.Q0 = 'Isotropic' by default.
%  - Type.Q is equal to one of 'Isotropic', 'Diag', 'Full'.  Type.Q = 'Isotropic' by default.
%  - Type.R is equal to one of 'Isotropic', 'Diag', 'Full'.  Type.R = 'Isotropic' by default.
% 'Model', followed by a struct Model of the parameters to initialize with.  Model = initialize_parameters(I, J) by default.
%
% ---Outputs---
% model:  learned MLDS parameters after running EM algorithm
% diagnostics:  struct of diagnostics such that
%  - diagnostics.log_likelihood is a vector of log-likelihoods for EM iterations
%  - diagnostics.model is a cell array of parameter values for EM iterations
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
I = size(X{1});
N = numel(X);
X = ten2vec(X);
M = numel(I);
if I(M) == 1
  M = M - 1;
  I = I(1:M);
end
J = set_optional_argument('J', I, varargin);
W = ten2vec(set_optional_argument('W', vec2ten(zeros(prod(I), N), I), varargin));
maxiter = set_optional_argument('MaxIter', 10, varargin);
epsilon = set_optional_argument('Epsilon', 1e-5, varargin);
Type = set_optional_argument('Type', struct('Q0', 'Isotropic', 'Q', 'Isotropic', 'R', 'Isotropic'), varargin);
model = ten2vec(set_optional_argument('Model', initialize_parameters(I, J), varargin));

LOGLI = true;
if (nargout < 2)
  LOGLI = false;
end

ratio = 1;
diff = inf;
iter = 0;
oldLogli = -inf;
while ((abs(diff) > epsilon) && (iter < maxiter) && (~(isTiny(model.Q0) || isTiny(model.Q) || isTiny(model.R))))
  oldmodel = model;
  iter = iter + 1;
  elapsed_time = tic;
  if (LOGLI)
    [mu, V, P, logli] = forward(X, model);
  else
    [mu, V, P] = forward(X, model);
  end
  [Ez, Ezz, Ez1z] = backward(mu, V, P, model);
  model = MLE_mlds(X, Ez, Ezz, Ez1z, model, Type);
  for n = 1:N
    xHAT = model.C*Ez{n};
    X(find(W(:,n)), n) = xHAT(find(W(:,n)));
  end

  if (LOGLI)
    logli = real(logli);
    diff = (logli - oldLogli);
    if (logli < oldLogli)
      warning('Log-likelihood decreases!');
    end
    ratio = abs(diff/logli);
    oldLogli = logli;
    fprintf('iteration = %d, logli = %d, time = %d\n', iter, logli, toc(elapsed_time));
  else
    fprintf('iteration = %d\n', iter);
  end

  diagnostics.log_likelihood(iter) = logli;
  diagnostics.model{iter} = vec2ten(model);
end
model = vec2ten(oldmodel);
end

function [t] = isTiny(sigma)
% test whether the matrix sigma is close to zero
  t = (norm(sigma, 1) < eps) || (any(diag(sigma) < eps));
end 
