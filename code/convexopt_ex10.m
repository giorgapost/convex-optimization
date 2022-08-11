% Author: Georgios Apostolakis
close all;
clear;
clc;

RANDOM = 1; %Constants
EXACT = 2;
MAT_CONSTR_METHODS = [RANDOM, EXACT];

rng(1);
%% Initialization
m = 5;   n = 10; %Dimensions of the problem
A = rand(m,n); x0 = ones(n,1)/n; x_true = rand(n,1);

for method = 1:length(MAT_CONSTR_METHODS)
    if MAT_CONSTR_METHODS(method) == RANDOM
        b{method}=rand(m,1);
    elseif MAT_CONSTR_METHODS(method) == EXACT
        b{method} = A*x_true + randn([m,1]);
    end
end

f = @(x, method) norm(A*x-b{method},1);
sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book
subgrad_f = @(x, method) A'*sgn(A*x-b{method});

%% Question 10.a: Solve with CVX
for method = 1:length(MAT_CONSTR_METHODS)
    cvx_begin quiet
        variable opt_pnt(n)
        minimize f(opt_pnt, method)
        subject to
            opt_pnt>=0;
            sum(opt_pnt)==1;
    cvx_end
    x_opt{method} = opt_pnt;
    f_opt(method) = cvx_optval;
end

%% Question 10.b: Subgradient algorithm with adaptive stepsize
for method = 1:length(MAT_CONSTR_METHODS)
    ff = @(x) f(x,method);
    subgrad_ff = @(x) subgrad_f(x,method);
    c = 1.001; %Loop termination parameter
    
    %Projection onto the unit simplex
    f_proj = @(x, mu) ones(n,1)' * max((x-mu*ones(n,1)), zeros(n,1)) - 1;
    mu_star = @(x) ex10_bisection(f_proj, x, min(x)-1, max(x)+1, 10^-10);
    proj = @(x) max(x-mu_star(x)*ones(n,1), zeros(n,1)); 
    
    theta = 1; L_f = norm(A,2)*sqrt(m);
    adaptive_step = @(x, k) (norm(subgrad_ff(x))~=0)*(sqrt(2*theta)) / (norm(subgrad_ff(x))*sqrt(k+1)) + (norm(subgrad_ff(x))==0)*(sqrt(2*theta)) / (L_f*sqrt(k+1)); %Adaptive step
    [x_subgr_alg{method}, f_subgr_alg{method}] = ex10_proj_subgrad_alg(x0, ff, subgrad_ff, adaptive_step, proj, f_opt(method), c);
end

%% Question 10.c: Mirror descent algorithm
for method = 1:length(MAT_CONSTR_METHODS)
    ff = @(x) f(x,method);
    subgrad_ff = @(x) subgrad_f(x,method);
    c = 1.001; %Loop termination parameter

	L_f = norm(abs(A')*ones(m,1), Inf);
    adaptive_step = @(x, k) (norm(subgrad_ff(x))~=0)*(sqrt(2)) / (norm(subgrad_ff(x), Inf)*sqrt(k+1)) + (norm(subgrad_ff(x))==0)*(sqrt(2)) / (L_f*sqrt(k+1)); %Adaptive step
    [x_mirror{method}, f_mirror{method}] = ex10_mirror_descent(x0, ff, subgrad_ff, adaptive_step, f_opt(method), c);
end

%% Question 10.d: Plots
for method = 1:length(MAT_CONSTR_METHODS)
    figure;
    semilogy(0:length(f_subgr_alg{method})-1, f_subgr_alg{method}-f_opt(method), 'DisplayName', 'Proj. subgrad.: f(x_k)-f_{opt}'); hold on;
    semilogy(0:length(f_mirror{method})-1, f_mirror{method}-f_opt(method), 'DisplayName', 'Mirror desc.: f(x_k)-f_{opt}');
    legend(); grid on; xlabel('Iterations (k)'); ylabel('f(x_k)-f_{opt}');
end