% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
rng(0);
m = 7;   n = 6; %Dimensions of the problem
A = rand(m,n); x_true = rand(n,1); b = A*x_true; 
d = rand(n,1); v1 = -100*rand(n,1); v2 = 100*rand(n,1); assert(sum(v1<=v2)==n);

f = @(x) 0.5*(x-d)'*(x-d);

%% Question 12.a
cvx_begin quiet  %Check that S1 is a nonempty set
    variable xx(n)
    minimize 0
    subject to
        A*xx<=b;
cvx_end
if ~strcmp(cvx_status, 'Solved')
    error('The set S1 is empty');
end

cvx_begin quiet %Obtain the optimal value f_opt
    variable x_cvx(n)
    minimize f(x_cvx)
    subject to
        A*x_cvx<=b;
cvx_end
f1_opt = cvx_optval;

prox_func = @(z) min(z,b);  %Solve the optimization problem with DPG
argmax_sol_func = @(y) d + A'*y;
[x_S1, f_S1] = ex12_DPG(A, f, 1, argmax_sol_func, prox_func, f1_opt, 10^-9);

figure;  %Plots
semilogy(0:length(f_S1)-1, abs(f_S1-f1_opt), 'DisplayName', 'DPG: f(x_k)-f_{opt}');
legend(); grid on; xlabel('Iterations (k)'); ylabel('f(x_k)-f_{opt}'); title('Q. 12.a');

%% Question 12.b
cvx_begin quiet  %Check that S2 is a nonempty set
    variable xx(n)
    minimize 0
    subject to
        A*xx == b;
        v1 <= xx;
        xx <= v2;
cvx_end
if ~strcmp(cvx_status, 'Solved')
    error('The set S2 is empty');
end

cvx_begin quiet %Obtain the optimal value f_opt
    variable x_cvx(n)
    minimize f(x_cvx)
    subject to
        A*x_cvx == b;
        v1 <= x_cvx;
        x_cvx <= v2;
cvx_end
f2_opt = cvx_optval;

A0 = [A;eye(n);eye(n)];  %Solve the optimization problem with DPG
prox_func = @(z) [b; max(z(m+1:m+n),v1); min(z(m+n+1:m+2*n),v2)];
argmax_sol_func = @(y) d + A0'*y;
[x_S2, f_S2] = ex12_DPG(A0, f, 1, argmax_sol_func, prox_func, f2_opt, 10^-9);

figure;  %Plots
semilogy(0:length(f_S2)-1, abs(f_S2-f2_opt), 'DisplayName', 'DPG: f(x_k)-f_{opt}');
legend(); grid on; xlabel('Iterations (k)'); ylabel('f(x_k)-f_{opt}'); title('Q. 12.b');

%% Question 12.c
cvx_begin quiet  %Check that S3 is a nonempty set
    variable xx(n)
    minimize 0
    subject to
        A*xx <= b;
        v1 <= xx;
        xx <= v2;
cvx_end
if ~strcmp(cvx_status, 'Solved')
    error('The set S3 is empty');
end

cvx_begin quiet %Obtain the optimal value f_opt
    variable x_cvx(n)
    minimize f(x_cvx)
    subject to
        A*x_cvx <= b;
        v1 <= x_cvx;
        x_cvx <= v2;
cvx_end
f3_opt = cvx_optval;

A0 = [A;eye(n);eye(n)];  %Solve the optimization problem with DPG
prox_func = @(z) [min(z(1:m),b); max(z(m+1:m+n),v1); min(z(m+n+1:m+2*n),v2)];
argmax_sol_func = @(y) d + A0'*y;
[x_S3, f_S3] = ex12_DPG(A0, f, 1, argmax_sol_func, prox_func, f3_opt, 10^-9);

figure;  %Plots
semilogy(0:length(f_S3)-1, abs(f_S3-f3_opt), 'DisplayName', 'DPG: f(x_k)-f_{opt}');
legend(); grid on; xlabel('Iterations (k)'); ylabel('f(x_k)-f_{opt}'); title('Q. 12.c');
