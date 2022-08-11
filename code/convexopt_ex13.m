% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
PIECEWISE = 1; 
SMOOTH = 2;
SIGNAL_SHAPES = [PIECEWISE, SMOOTH];

rng(0);
n = 100;
d_orig = zeros(n,length(SIGNAL_SHAPES));

for s=1:length(SIGNAL_SHAPES)
    if SIGNAL_SHAPES(s)==PIECEWISE
        pieces_n = [14, 8, 20, 18, 7, 9, 14, 10]; assert(sum(pieces_n)==n); %they must sum to n
        for i=1:length(pieces_n)
            d_orig(sum(pieces_n(1:i-1))+1:sum(pieces_n(1:i-1))+pieces_n(i), s) = randi([0,100]);
        end
    elseif SIGNAL_SHAPES(s)==SMOOTH
        tmp = 1:n;
        d_orig(:,s) = 10*sin(linspace(0,2*pi, n));
    end
end

d = d_orig + normrnd(0, 3, [n,length(SIGNAL_SHAPES)]);
f =  @(x,shape) 0.5*(x-d(:,shape))'*(x-d(:,shape));
F1 = @(x,lambda,shape) f(x,shape) + lambda*norm(x(1:end-1)-x(2:end),1);
F2 = @(x,lambda,shape) f(x,shape) + lambda/2*(x(1:end-1)-x(2:end))'*(x(1:end-1)-x(2:end));

%% Question 13.a - Solve with CVX various problems
LAMBDA_VALS = [0.1, 1, 10];

x_cvx = zeros(n,2,length(LAMBDA_VALS),length(SIGNAL_SHAPES));
f_cvx_l1 = zeros(length(LAMBDA_VALS),length(SIGNAL_SHAPES));

for s = 1:length(SIGNAL_SHAPES)
    for lam_idx = 1:length(LAMBDA_VALS)
        cvx_begin quiet
            variable x_opt(n)
            minimize F1(x_opt, LAMBDA_VALS(lam_idx), SIGNAL_SHAPES(s))
        cvx_end
        x_cvx(:,1,lam_idx,s) = x_opt;
        f_cvx_l1(lam_idx,s) = f(x_opt,SIGNAL_SHAPES(s));
        
        cvx_begin quiet
            variable x_opt(n)
            minimize F2(x_opt, LAMBDA_VALS(lam_idx), SIGNAL_SHAPES(s))
        cvx_end
        x_cvx(:,2,lam_idx,s) = x_opt;
    end
end

%% Question 13.a - Plots
for s = 1:length(SIGNAL_SHAPES)  
    figure;  hold on; %Plot F1 for various lambda
    plot(1:n, d_orig(:,s), 'DisplayName', 'Original signal');
    for lam_idx = 1:length(LAMBDA_VALS)
        plot(1:n, x_cvx(:,1,lam_idx,s), 'DisplayName', sprintf('Solution for lambda = %.1f',LAMBDA_VALS(lam_idx)));
    end
    legend; grid on; title('L1-Regularization'); xlabel('Dimension (i)'); ylabel('x_i');
    
    figure;  hold on;%Plot F2 for various lambda
    plot(1:n, d_orig(:,s), 'DisplayName', 'Original signal');
    for lam_idx = 1:length(LAMBDA_VALS)
        plot(1:n, x_cvx(:,2,lam_idx,s), 'DisplayName', sprintf('Solution for lambda = %.1f',LAMBDA_VALS(lam_idx)));
    end
    legend; grid on; title('L2-Regularization'); xlabel('Dimension (i)'); ylabel('x_i');
end

%% Question 13.b - Solve with DPG, FDPG the L1 problem
x_DPG = cell(length(LAMBDA_VALS), length(SIGNAL_SHAPES));
f_DPG = cell(length(LAMBDA_VALS), length(SIGNAL_SHAPES));
x_FDPG = cell(length(LAMBDA_VALS), length(SIGNAL_SHAPES));
f_FDPG = cell(length(LAMBDA_VALS), length(SIGNAL_SHAPES));

DD=eye(n)-circshift(eye(n),-1); DD=DD(1:n-1, :); normD=norm(DD,2);
D = @(x) x(1:end-1) - x(2:end);
D_T = @(y) [y;0] - [0;y];
L = norm(DD,2)^2/1;

for s = 1:length(SIGNAL_SHAPES)
    argmax_sol_func = @(y) d(:,s) + D_T(y);
    ff = @(x) f(x,SIGNAL_SHAPES(s));
    for lam_idx = 1:length(LAMBDA_VALS)
        prox_func = @(z) ex06_soft_thresh(z,LAMBDA_VALS(lam_idx)*L);
        [x_DPG{lam_idx,s}, f_DPG{lam_idx,s}] = ex13_DPG(n, ff, L, D, argmax_sol_func, prox_func, f_cvx_l1(lam_idx,s), 10^-3);
        [x_FDPG{lam_idx,s}, f_FDPG{lam_idx,s}] = ex13_FDPG(n, ff, L, D, argmax_sol_func, prox_func, f_cvx_l1(lam_idx,s), 10^-3);
    end
end

for s = 1:length(SIGNAL_SHAPES)  
    figure; %Plot for various lambda
    for lam_idx = 1:length(LAMBDA_VALS)
        semilogy(0:length(f_DPG{lam_idx,s})-1, abs(f_DPG{lam_idx,s}-f_cvx_l1(lam_idx,s)), 'DisplayName', sprintf('lambda=%.1f',LAMBDA_VALS(lam_idx)));
        hold on;
    end
    legend; grid on; ylabel('f(x_k)-f_{opt}'); xlabel('Iterations (k)'); 
    if SIGNAL_SHAPES(s)==PIECEWISE
        title('DPG - Piecewise signal');
    else
        title('DPG - Sinusoidal signal');
    end
end

for s = 1:length(SIGNAL_SHAPES)  
    figure; %Plot for various lambda
    for lam_idx = 1:length(LAMBDA_VALS)
        semilogy(0:length(f_FDPG{lam_idx,s})-1, abs(f_FDPG{lam_idx,s}-f_cvx_l1(lam_idx,s)), 'DisplayName', sprintf('lambda=%.1f',LAMBDA_VALS(lam_idx)));
        hold on;
    end
    legend; grid on; ylabel('f(x_k)-f_{opt}'); xlabel('Iterations (k)'); 
    if SIGNAL_SHAPES(s)==PIECEWISE
        title('FDPG - Piecewise signal');
    else
        title('FDPG - Sinusoidal signal');
    end
end