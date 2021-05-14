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
function D = err(X, Ntrain, model)
%
% compute mean-squared error.
%
% ---Inputs---
% X:  [N x 1] cell array representing a tensor time series
% Ntrain:  Postive integer indicating the number of epochs to train on.  Train on the first Ntrain tensors and test on the rest.
% model:  MLDS parameters
%
% ---Outputs---
% D:  [Ntest x 1] vector such that D(n) is the relative Euclidean error ||X_n - Xest_n|| / ||X_n||
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
N = numel(X);
I = size(X{1});
J = size(model.mu0);
p = prod(I);
M = numel(I);
if I(M) == 1
  M = M - 1;
  I = I(1:M);
end
J = J(1:M);
q = prod(J);
X = ten2vec(X);
Ntest = N - Ntrain;
model = ten2vec(model);
Xtrain = X(:,1:Ntrain);
Xtest = X(:,[1:Ntest]+Ntrain);

% estimate
[mu V P logli] = forward(Xtrain, model);
[Ez Ezz Ez1z] = backward(mu, V, P, model);
model.mu0 = model.A*Ez{Ntrain};
[ZEst ZZEst P logli] = forward(Xtest, model);
XEst = zeros(p,Ntest);
for i = 1:Ntest
  XEst(:,i) = model.C*ZEst{i};
end

% error
D = zeros(Ntest,1);
for i = 1:Ntest
  D(i) = norm(X(:,i+Ntrain)-XEst(:,i),'fro')/norm(X(:,i+Ntrain),'fro');
end
