close all;
clear;
clc;

m = 5;   n = 10; %Dimensions of the problem
rng(0);
for ptr = 1:2  %Repeated executions of the script
    %% Initialization
    A=rand(m,n);    b=rand(m,1);
    x0=zeros(n,1); x0 = x0 - A'*((A*A')\(A*x0-b)); %Projection onto Ax=b

    h = @(x) norm(x,1);
    sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book
    subgrad_h = @(x) sgn(x);
    grad_h_smooth = @(x,mu) (1/mu)*(x-ex06_soft_thresh(x, mu));

    %% Question 9.a: Solve with CVX
    cvx_begin quiet
        variable x_opt(n)
        minimize h(x_opt)
        subject to
            A*x_opt==b;
    cvx_end
    h_opt = cvx_optval;

    %% Questions 9.b, 9.c: Solve with S-FISTA, Projected subgradient alg.
    %S-FISTA algorithm
    epsilon = 10^-4;
    proj = @(x) x - A'*((A*A')\(A*x-b));
    [x_SFISTA, h_SFISTA] = ex09_SFISTA_alg(x0, h, grad_h_smooth, proj, h_opt, epsilon);

    %Subgradient algorithm with Polyak step
    polyak_step = @(x, k) (norm(subgrad_h(x))>0).*(h(x)-h_opt) / (norm( subgrad_h(x) )^2) + (norm(subgrad_h(x))==0)*1; %Polyak step
    [x_pol, h_pol] = ex09_subgrad_alg(x0, h, subgrad_h, polyak_step, proj, length(h_SFISTA));

    %% Question 9.d: Plots
    figure;
    semilogy(0:length(h_pol)-1, h_pol-h_opt, 'DisplayName', 'Proj. subgrad.: h(x_k)-h_{opt}'); hold on;
    semilogy(0:length(h_SFISTA)-1, h_SFISTA-h_opt, 'DisplayName', 'S-FISTA:  h(x_k)-h_{opt}');
    legend(); grid on; xlabel('Iterations (k)'); ylabel('h(x_k)-h_{opt}');
end