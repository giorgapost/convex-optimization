% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Definition of constants & functions, Setup of the problem

DIMENSIONS = [ 25,  5;  %(m, n) dimensions of the matrices
               25, 25;
                5, 25];
rng(0);
for dim_ptr = 1:size(DIMENSIONS,1)
    m = DIMENSIONS(dim_ptr,1);    n = DIMENSIONS(dim_ptr,2);

%%  Initialization
    A = rand(m,n);    b = rand(m,1);    x0 = rand(n,1);
    lambda = 1;
    
    sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book
    f = @(x) norm(A*x-b, 1);     subgrad_f = @(x) A'*sgn(A*x-b);
    grad_f_smooth = @(x, mu) (1/mu)*A'*(A*x-b-ex06_soft_thresh(A*x-b,mu));
    g = @(x) lambda*norm(x,1);   subgrad_g = @(x) lambda*sgn(x);
    F = @(x) f(x) + g(x);  subgrad_F = @(x) subgrad_f(x) + subgrad_g(x);
    grad_F_smooth = @(x,mu) grad_f_smooth(x,mu) + (lambda/mu)*(x-ex06_soft_thresh(x,mu));
    
    %Solve with CVX
    cvx_begin quiet
        variable x_opt(n)
        minimize F(x_opt)
    cvx_end
    F_opt = cvx_optval;
    
%%  Projected subgradient algorithm
    c = 1.01;  %Loop termination criterion
    
    %Projected subgradient algorithm with Polyak step
    polyak_step = @(x, k) (norm(subgrad_F(x))>0).*(F(x)-F_opt) / (norm( subgrad_F(x) )^2) + (norm(subgrad_F(x))==0)*1;
    [x_pol, F_pol] = ex08_proj_subgrad_alg(x0, F, subgrad_F, polyak_step, F_opt, c);

    %Projected subgradient algorithm with Dynamic step
    L_F = norm(A,2)*sqrt(m)+lambda*sqrt(n);
    dynamic_step = @(x, k) (norm(subgrad_F(x))>0).*1/(norm(subgrad_F(x))*sqrt(k+1)) + (norm(subgrad_F(x))==0)*(1/L_F);
    [x_dyn, F_dyn] = ex08_proj_subgrad_alg(x0, F, subgrad_F, dynamic_step, F_opt, c);
    
%%  Proximal subgradient algorithm with Dynamic step
    [x_prox, F_prox] = ex08_prox_subgrad_alg(x0, F, subgrad_f, dynamic_step, lambda, F_opt, c);

%%  S-FISTA algorithm - Smoothening only f
    epsilon = (c-1)*F_opt;
    alpha = norm(A,2)^2;
    beta = m/2;
    [x_SFISTA_c, F_SFISTA_c] = ex08_SFISTA_alg(x0, F, grad_f_smooth, 0, alpha, beta, lambda, F_opt, epsilon);
    
%%  S-FISTA algorithm - Smoothening both f and g
    alpha = norm(A,2)^2+lambda;
    beta = m/2 + lambda*n/2;
    [x_SFISTA_d, F_SFISTA_d] = ex08_SFISTA_noprox_alg(x0, F, grad_F_smooth, 0, alpha, beta, F_opt, epsilon);

%%  Plots
    F_pol_best = zeros(size(F_pol));
    for k=1:length(F_pol)
        F_pol_best(k) = min(F_pol(1:k));
    end
    
    F_dyn_best = zeros(size(F_dyn));
    for k=1:length(F_dyn)
        F_dyn_best(k) = min(F_dyn(1:k));
    end
    
    F_prox_best = zeros(size(F_prox));
    for k=1:length(F_prox)
        F_prox_best(k) = min(F_prox(1:k));
    end

    figure;
    semilogy(0:length(F_pol)-1, F_pol-F_opt, 'DisplayName', 'Polyak: F(x_k)-F_{opt}'); hold on;
    semilogy(0:length(F_dyn)-1, F_dyn-F_opt, 'DisplayName', 'Dynamic: F(x_k)-F_{opt}');
    semilogy(0:length(F_pol_best)-1, F_pol_best-F_opt, 'DisplayName', 'Polyak: F^{best}_k-F_{opt}');
    semilogy(0:length(F_dyn_best)-1, F_dyn_best-F_opt, 'DisplayName', 'Dynamic: F^{best}_k-F_{opt}');    
    title(sprintf('m=%d, n=%d',m,n)); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');
    
    figure;
    semilogy(0:length(F_pol)-1, F_pol-F_opt, 'DisplayName', 'Polyak: F(x_k)-F_{opt}'); hold on;
    semilogy(0:length(F_dyn)-1, F_dyn-F_opt, 'DisplayName', 'Dynamic: F(x_k)-F_{opt}');
    semilogy(0:length(F_prox)-1, F_prox-F_opt, 'DisplayName', 'Proximal subgr.: F(x_k)-F_{opt}');
    semilogy(0:length(F_pol_best)-1, F_pol_best-F_opt, 'DisplayName', 'Polyak: F^{best}_k-F_{opt}');
    semilogy(0:length(F_dyn_best)-1, F_dyn_best-F_opt, 'DisplayName', 'Dynamic: F^{best}_k-F_{opt}');    
    semilogy(0:length(F_prox_best)-1, F_prox_best-F_opt, 'DisplayName', 'Proximal subgr.: F^{best}_k-F_{opt}');    
    title(sprintf('m=%d, n=%d',m,n)); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');

    figure;
    semilogy(0:length(F_SFISTA_c)-1, F_SFISTA_c-F_opt, 'DisplayName', 'Only f is smoothened');
    title(sprintf('S-FISTA, m=%d, n=%d',m,n)); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');
    
    figure;
    semilogy(0:length(F_SFISTA_c)-1, F_SFISTA_c-F_opt, 'DisplayName', 'Only f is smoothened'); hold on;
    semilogy(0:length(F_SFISTA_d)-1, F_SFISTA_d-F_opt, 'DisplayName', 'Both f and g are smoothened');
    title(sprintf('S-FISTA, m=%d, n=%d',m,n)); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');

end