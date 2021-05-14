% This script fits both the MLDS and the LDS to synthetic MLDS data 
% and compares the prediction accuracies and log-likelihoods.

% add the path to dynammo
addpath('dynammo');

% data
I = [7 11];
J = [3 5];
Type = struct('Q0','Isotropic','Q','Isotropic','R','Isotropic');
model = ten2vec(initialize_parameters(I, J));
M = numel(I);
N = 1100;
Ntrain = 1000;
Ntest = N - Ntrain;
Z = zeros(prod(J), N);
X = zeros(prod(I), N);
for n = 1:N
  if n == 1
    Z(:,n) = mvnrnd(model.mu0, model.Q0);
  else
    Z(:,n) = mvnrnd(model.A * Z(:,n-1), model.Q);
  end
  X(:,n) = mvnrnd(model.C * Z(:,n), model.R);
end
Z = vec2ten(Z, J);
X = vec2ten(X, I);
X_vectorized = vec2ten(ten2vec(X), prod(I));
model = vec2ten(model);

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

% true
X = ten2vec(X);
[mu_true V_true P_true l_true] = forward(X(:,1:Ntrain), ten2vec(model));
X = vec2ten(X, I);
clear A C mu_true V_true P_true
err_true = err(X, Ntrain, model);

disp(' ')

% plot results
disp('Plotting results...')
f = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
hold on;
T = [1:Ntest]+Ntrain; % time ticks
plot(T, err_lds, 'Color', 'blue');
plot(T, err_mlds, 'Color', 'red');
plot(T, err_true, 'Color', [0 .75 0], 'LineStyle', '--', 'LineWidth', 2);
hold off;
legend('LDS', 'MLDS', 'true');
xlim([1 Ntest] + Ntrain);
xlabel('Time slice');
ylabel('Error');

subplot(1,2,2);
numiter = max(numel(diagnostics_lds.log_likelihood), ...
              numel(diagnostics_mlds.log_likelihood));
hold on;
plot(diagnostics_lds.log_likelihood,'Color','blue');
plot(diagnostics_mlds.log_likelihood,'Color','red');
plot(l_true*ones(numiter,1), 'Color', [0 .75 0], 'LineStyle', '--', 'LineWidth', 2);
hold off;
legend('LDS','MLDS','true');
xlim([1,numiter]);
xlabel('Number of EM iterations');
ylabel('Log-likelihood');
