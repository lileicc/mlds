% This source code is (c) Copyright by Lei Li, Mark Rogers.
% All rights preserved.
%
% Permission is granted to use it for non-profit purposes,
% including research and teaching. For-profit use requires
% the express consent of the author (leili@cs.berkeley.edu).
%
% Details in the following paper:
%    Mark Rogers, Lei Li and Stuart J. Russell (2013),
%      "Multilinear Dynamical Systems for Tensor Time Series",    
%      In Advances in Neural Information Processing Systems 26. 
%
function [modelNEW] = MLE_mlds(X, Ez, Ezz, Ez1z, model, Type)
%
% Perform M-step of the EM algorithm
%
% ---Inputs---
% X:  prod(I) x N matrix, where X(:,n) is the vectorization of the nth tensor of the tensor time series
% Ez:  cell array such that Ez{n} is the conditional expectation of vec(Z_n) given (X_1, ..., X_N)
% Ezz:  cellarray such that Ezz{n} is the conditional expectation of vec(Z_n)*vec(Z_n)' given (X_1, ..., X_N)
% Ez1z:  cell array such that Ez1z{n} is the conditional expectation of vec(Z_n+1)*vec(Z_n)' given (X_1, ..., X_N)
% model:  "old" struct of the vectorized MLDS parameters
% Type.Q0:  either 'Isotropic', 'Diag', or 'Full'
% Type.Q:  either 'Isotropic', 'Diag', or 'Full'
% Type.R:  either 'Isotropic', 'Diag', or 'Full'
%
% ---Outputs---
% modelNEW: "new" struct of the vectorized MLDS parameters
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
if strcmp(class(model.cellC),'cell')
  M = numel(model.cellC);
  I = zeros(1,M);
  J = zeros(1,M);
  for m = 1:M
    I(m) = size(model.cellC{m},1);
    J(m) = size(model.cellC{m},2);
  end
else
  M = 1;
  I = size(model.cellC,1);
  J = size(model.cellC,2);
end
N = size(X, 2);
P = size(X, 1);
H = size(Ez{1}, 1);
Sz1z = zeros(H, H);
Szz = zeros(H, H);
Sxz = zeros(P, H);
for n = 1:(N-1)
  Sz1z = Sz1z + Ez1z{n};
end
for n = 1:N
  Szz = Szz + Ezz{n};
  Sxz = Sxz + X(:,n) * Ez{n}';
end
SzzN = Szz - Ezz{N};

% mu0, Q0
model.mu0 = Ez{1};
model.Q0 = Ezz{1} - Ez{1} * Ez{1}';
switch Type.Q0
case 'Diag'
  model.Q0 = diag(diag(model.Q0));
case 'Full'
  model.Q0 = Ezz{1} - Ez{1}*Ez{1}';
case 'Isotropic'
  model.Q0 = diag(repmat(trace(model.Q0) / H, H, 1));
end

% A
A_unfactorized = Sz1z / SzzN;
if M == 1
  model.cellA = A_unfactorized;
else
  switch Type.Q
  case 'Isotropic'
    OMEGA = model.Q(1)^(-1);
  case 'Diag'
    OMEGA = diag(model.Q) .^ (-1);
  case 'Full'
    OMEGA = model.Q \ eye(H);
  end
  PSI = zeros(H);
  PHI = zeros(H);
  for n = 1:N-1
    PSI = PSI + Ezz{n};
    PHI = PHI + Ez1z{n};
  end
  model.cellA = update_multilinear_operator(model.cellA, OMEGA, PSI, PHI, Type.Q);
end
model.A = mkron(model.cellA);

% Q
switch Type.Q
case 'Diag'
  model.Q = diag((diag(Szz) - diag(Ezz{1}) - 2 * diag(model.A * Sz1z') + diag(model.A * SzzN * model.A')) / (N-1));
case 'Full'
  tmp = model.A * Sz1z';
  model.Q = (Szz - Ezz{1} - tmp - tmp' + model.A * SzzN * model.A') / (N-1);
case 'Isotropic'
  delta = (trace(Szz) - trace(Ezz{1}) - 2 * trace(model.A * Sz1z') + trace(model.A * SzzN * model.A')) / (N-1) / H;
  model.Q = diag(repmat(delta, H, 1));
end

% C
C_unfactorized = Sxz / Szz;
if M == 1
  model.cellC = C_unfactorized;
else
  switch Type.R
  case 'Full'
    OMEGA = model.R \ eye(P);
  case 'Diag'
    OMEGA = diag(model.R) .^ (-1);
  case 'Isotropic'
    OMEGA = model.R(1)^(-1);
  end
  PSI = zeros(H);
  PHI = zeros(P,H);
  for n = 1:N
    PSI = PSI + Ezz{n};
    PHI = PHI + X(:,n)*Ez{n}';
  end
  model.cellC = update_multilinear_operator(model.cellC, OMEGA, PSI, PHI, Type.R);
end
model.C = mkron(model.cellC);

% R
switch Type.R 
case 'Diag'
  model.R = diag((diag(X * X') - 2 * diag(model.C * Sxz') + diag(model.C * Szz * model.C')) / N);
case 'Full'
  tmp = model.C * Sxz';
  model.R = (X * X' - tmp - tmp' + model.C * Szz * model.C') / N;
case 'Isotropic'
  delta = (trace(X * X') - 2 * trace(model.C * Sxz') + trace(model.C * Szz * model.C')) / N / P;
  model.R = diag(repmat(delta, P, 1));
end
  
modelNEW = model;  
