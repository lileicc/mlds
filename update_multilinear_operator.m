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
function BNEW = update_multilinear_operator(B, OMEGA, PSI, PHI, covariance_type)
%
% Perform M-step to update the multilinear operator B during the EM algorithm. B can either be the transition multilinear operator A or the projection multilinear operator C of MLDS.
%
% ---Inputs---
% B:  "old" value of the multilinear operator B, which is represented by a cell array of projection matrices
% OMEGA, PSI, PHI:  matrices constant with respect to B.  See /mlds/MLE_mlds.m.
% covariance_type:  either 'Isotropic', 'Diag', or 'Full'
%
% ---Outputs---
% BNEW: "new" value of the multilinear operator B, which is represented by a cell array of projection matrices
%
% @author: Mark Rogers (markrogersjr@berkeley.edu)
% @last modified date: 2013/12/13
%
  M = numel(B);
  I = zeros(1,M);
  J = zeros(1,M);
  for m = 1:M
    I(m) = size(B{m},1);
    J(m) = size(B{m},2);
  end
  sizeB = size(B);
  if ~strcmp(covariance_type,'Diag')
    OMEGA = (OMEGA + OMEGA')/2;
  end
  PSI = (PSI + PSI')/2;
  f = @(z) f(z, OMEGA, PSI, PHI, covariance_type, I, J);
  df = @(z) df(z, OMEGA, PSI, PHI, covariance_type, I, J);
  BNEW = tensor(descend(f,df,0,vector(B)),I,J);
end

%---------------------------------------------------------
function y = f(z, OMEGA, PSI, PHI, covariance_type, I, J)
  B = mkron(tensor(z, I,J));
  temp = PHI*(B');
  temp = B*PSI*B' - temp - temp';
  y = traceprod(OMEGA, temp, covariance_type);
end

function g = df(z, OMEGA, PSI, PHI, covariance_type, I, J)
  M = numel(I);
  B = tensor(z, I, J);
  pIJ = I .* J;
  H = PSI*mkron(B)'-PHI';
  if strcmp(covariance_type,'Diag')
    H = H*diag(OMEGA);
  else
    H = H*OMEGA;
  end
  H = H';
  g = zeros(sum(pIJ),1);
  for m = 1:M
    notm = find([1:M] ~= m);
    Bnotm = mkron(subcell(B,notm));
    Hm = mat2ten(H,[I J]);
    Hm = ten2mat(permute(Hm,[notm,m,notm+M,m+M]));
    G = zeros([I(m) J(m)]);
    for i = 1:I(m)
      for j = 1:J(m)
        G(i,j) = sum(vec(Bnotm .* Hm(...
          [1:prod(I(notm))] + prod(I(notm))*(i-1),...
          [1:prod(J(notm))] + prod(J(notm))*(j-1))));
      end
    end
    g([1:pIJ(m)] + sum(pIJ(1:m-1))) = 2*G(:);
  end
end

%-------------------------------------------------------------------
function v = vector(T)                                                   
  M = numel(T);
  sizes = zeros(M,1);
  for m = 1:M
    sizes(m) = numel(T{m});
  end
  v = zeros(sum(sizes),1);
  for m = 1:M
    v([1:sizes(m)]+sum(sizes(1:(m-1)))) = T{m}(:);
  end
end

%-------------------------------------------------------------------
function T = tensor(v,I,J)
  M = numel(I);
  sizes = zeros(M,1);
  for m = 1:M
    sizes(m) = I(m)*J(m);
  end
  T = cell(M,1);
  for m = 1:M
    T{m} = reshape(v([1:sizes(m)]+sum(sizes(1:(m-1)))),[I(m) J(m)]);
  end
end
