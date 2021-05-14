function [model, LL] = learn_lds(X, varargin)
% learning model parameters for Linear Dynamical Systems (LDS), also known
% as Kalman Filters.
% Linear Dynamical Systems are described by the following equations:
% z_1 = mu0 + w_1
% z_n = A z_{n-1} + w_n
% x_n = C x_n + v_n
% w_1 ~ N(0, Q0)
% w_n ~ N(0, Q)
% v_n ~ N(0, R)
%
% Args:
%   X: M * N matrix, M is number of sequences, N is the time duration.
%
% Optional Args:
%   'Hidden', followed by an integer indicating the number of hidden
%   variables.
%   'MaxIter', followed by an integer indicating the number of max
%   iterations.
%   'Model', followed by a struct denoting a model to start with (see
%   below).
%   'Fixed', followed by a model with some fieds of (A, C, Q, R, mu0, Q0), 
%           do not learn, use the provided parameters 
%  covariance options:
%   'DiagQ0', if presented, will learn a diagonal covariance Q0
%   'DiagQ', if presented, will learn a diagonal covariance Q
%   'DiagR', if presented, will learn a diagonal covariance R
%   'IsotropicQ0', if presented, will learn a diagonal and isotropic Q0
%   'IsotropicQ', if presented, will learn a diagonal and isotropic Q
%   'IsotropicR', if presented, will learn a diagonal and isotropic R
%   'FullQ0', if presented, will learn any possible Q0
%   'FullQ', if presented, will learn any possible Q
%   'FullR', if presented, will learn any possible R
%  Note these options could not coexist for the same covariance matrix.
%  Default (no args given) the algorithm will learn with H=M, MaxIter=10,
%  diagonal and isotropic Q0, Q, R.
%
% Returns:
%   model: a struct with the following attributes:
%     A: transition matrix (also named as F in papers), H * H
%     C: transmission  matrix(also called projection matrix, named as G
% in papers), M * H
%     Q: transition covariance, H * H 
%     R: transmission covariance, M * M
%     mu0: initial states (also named as z_0 sometimes), H * 1
%     Q0: initial covariance, H * H
%
% Example:
% t = 1:100;
% x1 = sin(2 * pi * t / 50);
% x2 = sin(2 * pi * t / 50 + pi / 4);
% X = [x1; x2];
% model = learn_lds(X, 'Hidden', 2, 'MaxIter', 100);
%
% $Author$@cs.berkeley.edu
% $Date$
% $Rev$
%

N = size(X, 2);
M = size(X, 1);

% get number of hidden variables
a = find(strcmp('MaxIter', varargin), 1);
if (isempty(a))
  maxIter = 10;
else
  maxIter = varargin{a+1};
end

% get number of hidden variables
a = find(strcmp('Hidden', varargin), 1);
if (isempty(a))
  H = M;
else
  H = varargin{a+1};
end

% get the initial model
a = find(strcmp('Model', varargin), 1);
if (isempty(a))
  model.A = eye(H, H) + randn(H, H);
  model.C = eye(M, H) + randn(M, H);
  model.Q = eye(H, H);
  model.R = eye(M, M);
  model.mu0 = randn(H, 1);
  model.Q0 = model.Q;
else
  model = varargin{a+1};
end

% if there is fixed parameter for part of the model
a = find(strcmp('Fixed', varargin), 1);
if (~isempty(a))
  modelinit = varargin{a+1};
  if (isfield(modelinit, 'A'))
    model.A = modelinit.A;
  end
  if (isfield(modelinit, 'C'))
    model.C = modelinit.C;
  end
  if (isfield(modelinit, 'Q'))
    model.Q = modelinit.Q;
  end
  if (isfield(modelinit, 'R'))
    model.R = modelinit.R;
  end
  if (isfield(modelinit, 'mu0'))
    model.mu0 = modelinit.mu0;
  end
  if (isfield(modelinit, 'Q0'))
    model.Q0 = modelinit.Q0;
  end
end

% plot eigen value trajectory of A
a = find(strcmp('ShowEigen', varargin), 1);
if (isempty(a))
  SHOWEIGEN = false;
else
  SHOWEIGEN = true;
  plotEigen(model.A, '.');
  hold on;
end

LOGLI = true;
if (nargout < 2)
  LOGLI = false;
end

CONV_BOUND = 1e-5;

ratio = 1;
diff = 1;
iter = 0;
oldLogli = -inf;

oldmodel = model;
while ((ratio > CONV_BOUND || diff > CONV_BOUND) && (iter < maxIter) && (~ (isTiny(model.Q0) || isTiny(model.Q) || isTiny(model.R))))
  oldmodel = model;
  iter = iter + 1;
  if (LOGLI)
    [mu, V, P, logli] = forward(X, model, varargin{:});
  else 
    [mu, V, P] = forward(X, model, varargin{:});
  end
  [Ez, Ezz, Ez1z] = backward(mu, V, P, model);
  model = MLE_lds(X, Ez, Ezz, Ez1z, varargin{:});
  if (LOGLI)
    logli = real(logli);
    diff = (logli - oldLogli);
    if (logli < oldLogli)
      warning('Loglikelihood decreases!');
    end
    ratio = abs(diff/logli) ;
    LL(iter) = logli;
    oldLogli = logli;
    fprintf('iteration = %d, logli = %d\n', iter, logli);
  else
    fprintf('iteration = %d\n', iter);
  end
  if (SHOWEIGEN)
    plotEigen(model.A, '.');
    drawnow;
  end
end
model = oldmodel;

if (SHOWEIGEN)
  plotEigen(model.A, 'kx');
  drawnow;
end

function [t] = isTiny(sigma)
% test whether the matrix sigma is close to zero
%
t = (norm(sigma, 1) < eps) || (any(diag(sigma) < eps));