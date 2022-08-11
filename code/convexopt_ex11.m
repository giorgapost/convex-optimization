% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
rng(0);
m = 5;   n = 500; %Dimensions of the problem
A = rand(m,n); b = rand(m,1); c = rand(n,1);

f = @(x) c'*x;
F = @(x, z) f(x) + (min(z)<0)*realmax;

cvx_begin quiet
    variable x_star(n)
    minimize f(x_star)
    subject to
        x_star>=0;
        A*x_star==b;
cvx_end
F_opt = cvx_optval; %the opt. value of F is the same with that of f

%% ADMM
epsilon = 10^(-4); %Loop termination parameter
x0 = 10*rand(n,1); z0 = 10*rand(n,1); p=20;

upd_x = @(x,z,y) (z-(1/p)*(y+c))- (A'/(A*A'))*(A*(z-(1/p)*(y+c))) + (A'/(A*A'))*b;
upd_z = @(x,z,y) max(x+(1/p)*y, 0);
upd_y = @(x,z,y) y+p*(x-z);

[x_ADMM, ~, F_ADMM] = ex11_ADMM(x0, z0, upd_x, upd_z, upd_y, F, F_opt, epsilon);

figure;
semilogy(0:length(F_ADMM)-1, abs(F_ADMM-F_opt), 'DisplayName', 'ADMM: |f(x_k)-f_{opt}|');
legend(); grid on; xlabel('Iterations (k)'); ylabel('|f(x_k)-f_{opt}|');