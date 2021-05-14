% This script fits both the MLDS and the LDS to the SST dataset
% and compares the prediction accuracies and log-likelihoods.

% add the path to dynammo
addpath('dynammo');

% data
load data/sst
I = size(X{1});
J = [3 3];
Ntrain = 1800;
Type = struct('Q0','Isotropic','Q','Isotropic','R','Isotropic');
M = numel(I);
N = numel(X);
Ntest = N - Ntrain;
X_vectorized = vec2ten(ten2vec(X), prod(I));

% LDS
disp('Fitting LDS...')
sizes = zeros(prod(I), 1);
for i = 1:prod(I)
  sizes(i) = number_of_parameters(prod(I), i, Type);
end
clear i
J_lds = find(sizes >= number_of_parameters(I, J, Type),1);
clear i sizes
[model_lds diagnostics_lds] = learn_mlds(subcell(X_vectorized, 1:Ntrain), 'J', J_lds);
err_lds = err(X, Ntrain, model_lds);

disp(' ')

% MLDS
disp('Fitting MLDS with matching number of parameters...')
J_mlds = prod(J);
[model_mlds diagnostics_mlds] = learn_mlds(subcell(X, 1:Ntrain), 'J', J);
err_mlds = err(X, Ntrain, model_mlds);

disp(' ')

% plot results
disp('Plotting results...')
f = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
hold on;
T = [1:Ntest]+Ntrain; % time ticks
plot(T, err_lds, 'Color', 'blue');
plot(T, err_mlds, 'Color', 'red');
hold off;
legend('LDS', 'MLDS');
xlim([1 Ntest] + Ntrain);
xlabel('Time slice');
ylabel('Error');

subplot(1,2,2);
numiter = max(numel(diagnostics_lds.log_likelihood), ...
              numel(diagnostics_mlds.log_likelihood));
hold on;
plot(diagnostics_lds.log_likelihood,'Color','blue');
plot(diagnostics_mlds.log_likelihood,'Color','red');
hold off;
legend('LDS','MLDS');
xlim([1,numiter]);
xlabel('Number of EM iterations');
ylabel('Log-likelihood');
