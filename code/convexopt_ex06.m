% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
rng(0);
n = 128;
s = 4; %s << n
m = 40; %m >= 2s*log(n) << n

A = rand(m,n);
xs_original = zeros(n,1);
xs_permuted_idxs = randperm(n);
xs_original(xs_permuted_idxs(1:s)) = 10*rand(s,1);
b = A*xs_original;

f = @(x) 0.5*(A*x-b)'*(A*x-b);
grad_f = @(x) (A')*(A*x-b);
g = @(x, lambda) lambda*norm(x,1);
f1 = @(x,lambda) f(x) + g(x,lambda);
f2 = @(x,lambda) f(x) + lambda*x'*x;

lambda_vals = [0.01, 0.1, 1, 10];
x1_star = zeros(n,length(lambda_vals));
x2_star = zeros(n,length(lambda_vals));
f1_opt = zeros(1, length(lambda_vals));
f2_opt = zeros(1, length(lambda_vals));

%% Question 6.a
for l = 1:length(lambda_vals)
    cvx_begin quiet
        variable x1(n)
        minimize f1(x1, lambda_vals(l))
    cvx_end
    x1_star(:,l) = x1;
    f1_opt(l) = cvx_optval;
    
    cvx_begin quiet
        variable x2(n)
        minimize f2(x2, lambda_vals(l))
    cvx_end
    x2_star(:,l) = x2;
    f2_opt(l) = cvx_optval;
end
x_ex06a = [xs_original, x1_star, x2_star];

%% Question 6.b
x0 = rand(n,1);
c = 1.01;
Lf = max(eig(A'*A));
sigma = min(eig(A'*A));

for l = 2:length(lambda_vals) %we omit the case where lambda=0.01
    %ISTA algorithm
    [x_ISTA{l}, f_ISTA{l}] = ex06_ISTA_alg(x0, f1, grad_f, Lf, lambda_vals(l), f1_opt(l), c);
    fprintf('Completed ISTA\n');
    
    %FISTA algorithm
    [x_FISTA{l}, f_FISTA{l}] = ex06_FISTA_alg(x0, f1, grad_f, Lf, lambda_vals(l), f1_opt(l), c);
    fprintf('Completed FISTA\n');
end

%Plots
for l = 2:length(lambda_vals) %we omit the case where lambda=0.01
    figure;
    semilogy(0:length(f_ISTA{l})-1, f_ISTA{l}-f1_opt(l), 'DisplayName', 'ISTA'); hold on;
    semilogy(0:length(f_FISTA{l})-1, f_FISTA{l}-f1_opt(l), 'DisplayName', 'FISTA'); 
    title(sprintf('Lambda=%.2f', lambda_vals(l))); legend();
    grid on; xlabel('Iterations (k)'); ylabel('f_1(x_k)-f_{1(opt)}');
end
