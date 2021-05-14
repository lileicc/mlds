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
function x = descend(f, df, ddf, x0)
%
% standard gradient descent algorithm.  Computes local minimum of f nearest the point x0.
%
% ---Inputs---
% f:  function from R^n to R, where n is the dimensionality of x0
% df:  gradient function of f
% ddf:  Hessian of f
% x0:  vector in R^n. Initial point
%
% ---Outputs---
% x:  vector in R^n.  Location of minimum closest to x_0 with respect to f.
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
epsilon = 1e-18;
maxiter = 1e12;
learning_rate = .8;
step_size = learning_rate;

numiter = 0;
x = x0;
f_x = inf;
f_x0 = -inf;
while abs(f_x - f_x0)/step_size > epsilon
  elapsed_time = tic;
  numiter = numiter + 1;
  x0 = x;
  f_x0 = f(x0);
%  x = x0 - step_size * pinv(ddf(x0)) * df(x0);
  x = x0 - step_size * df(x0);
  f_x = f(x);
  while f_x > f_x0
    step_size = step_size * learning_rate;
%    x = x0 - step_size * pinv(ddf(x0)) * df(x0);
    x = x0 - step_size * df(x0);
    f_x = f(x);
  end
  %disp(['  iter=' num2str(numiter) ' df=' num2str((f(x)-f(x0))/step_size,'%d') ' step=' num2str(step_size,'%d') ' dt=' num2str(toc(elapsed_time),'%d')]);
end
