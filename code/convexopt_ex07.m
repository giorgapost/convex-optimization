% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
WEAK = 0.2; %std.deviation of weak noise level
STRONG = 2.5; %std. deviation of strong noise level
NOISE_LEVEL = [WEAK, STRONG];

rng(0);
m = 20;
n = 15;
pieces_n = [4, 3, 5, 3]; %they must sum to n
x_pwc = zeros(n,1);
for i=1:length(pieces_n)
    x_pwc(sum(pieces_n(1:i-1))+1:sum(pieces_n(1:i-1))+pieces_n(i)) = randi([0,10]);
end
A = rand(m,n);  b = zeros(m,2);  x0 = rand(n,1);
for lvl = 1:length(NOISE_LEVEL)
    b(:,lvl) = A*x_pwc + normrnd(0,NOISE_LEVEL(lvl),m,1);
end

f = @(x, noise_lvl) 0.5*(A*x-b(:,noise_lvl))'*(A*x-b(:,noise_lvl));
grad_f = @(x, noise_lvl) A'*(A*x-b(:,noise_lvl));
p = 1;
D = @(x) p*(x(1:n-1) - x(2:n));
D_T = @(y) [y;0] - [0;y];
DD=p*(eye(n)-circshift(eye(n),-1)); DD=DD(1:n-1, :); normD=norm(DD,2);
F = @(x, noise_lvl) f(x,noise_lvl) + norm(D(x), 1);
sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book
subgrad_F = @(x, noise_lvl) A'*(A*x-b(:,noise_lvl)) + D_T(sgn(D(x)));

%% Question 7.a
for lvl = 1:length(NOISE_LEVEL)
    cvx_begin quiet
        variable x_cvx(n)
        minimize f(x_cvx, lvl)
    cvx_end
    xf_opt(:,lvl) = x_cvx;
    f_opt(lvl) = cvx_optval;
end

%% Question 7.b
for lvl = 1:length(NOISE_LEVEL)
    cvx_begin quiet
        variable x_cvx(n)
        minimize F(x_cvx, lvl)
    cvx_end
    xF_opt(:,lvl) = x_cvx;
    F_opt(lvl) = cvx_optval;
end

x_ex07ab = [x_pwc, xf_opt, xF_opt];

%% Question 7.c, 7.d
for lvl = 1:length(NOISE_LEVEL)
    
    %7c. S-FISTA algorithm
    FF = @(x) F(x,lvl);
    subgrad_FF = @(x) subgrad_F(x,lvl);
    grad_ff = @(x) grad_f(x,lvl);
    Lf = norm(A'*A, 2);
    epsilon = 0.001;
    [x_SFISTA{lvl}, F_SFISTA{lvl}] = ex07_SFISTA_alg(x0, FF, grad_ff, Lf, D, D_T, normD, F_opt(lvl), epsilon);
    
    %7d. Subgradient algorithm with Polyak step
    polyak_step = @(x, k) (norm(subgrad_FF(x))>0).*(FF(x)-F_opt(lvl)) / (norm( subgrad_FF(x) )^2) + (norm(subgrad_FF(x))==0)*1; %Polyak step
    [x_pol{lvl}, F_pol{lvl}] = ex07_subgrad_alg(x0, FF, subgrad_FF, polyak_step, length(F_SFISTA{lvl}));
end

%% Question 7.e: Plots
for lvl = 1:length(NOISE_LEVEL)
    F_pol_best{lvl} = zeros(size(F_pol{lvl}));
    for k=1:length(F_pol{lvl})
        F_pol_best{lvl}(k) = min(F_pol{lvl}(1:k));
    end
    
    figure;
    semilogy(0:length(F_SFISTA{lvl})-1, F_SFISTA{lvl}-F_opt(lvl), 'DisplayName', 'S-FISTA'); hold on;
    title(sprintf('Noise Std. Dev. = %.2f', NOISE_LEVEL(lvl))); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');
    
    figure;
    semilogy(0:length(F_pol{lvl})-1, F_pol{lvl}-F_opt(lvl), 'DisplayName', 'F(x_k)-F_{opt}'); hold on;
    semilogy(0:length(F_pol_best{lvl})-1, F_pol_best{lvl}-F_opt(lvl), 'DisplayName', 'F^{best}_k-F_{opt}');
    title(sprintf('Subgrad. desc. - Polyak, Noise Std. Dev. = %.2f', NOISE_LEVEL(lvl))); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');
    
    figure;
    semilogy(0:length(F_SFISTA{lvl})-1, F_SFISTA{lvl}-F_opt(lvl), 'DisplayName', 'S-FISTA'); hold on;
    semilogy(0:length(F_pol{lvl})-1, F_pol{lvl}-F_opt(lvl), 'DisplayName', 'Subgrad. desc. - Polyak');
    title(sprintf('Noise Std. Dev. = %.2f', NOISE_LEVEL(lvl))); legend();
    grid on; xlabel('Iterations (k)'); ylabel('F(x_k)-F_{opt}');
end