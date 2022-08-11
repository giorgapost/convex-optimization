% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
rng(0);
n = 2;  %dimension of the problem
x = rand(n,1);

%% Computing the projection
f = @(mu) ones(n,1)' * max((x-mu*ones(n,1)), zeros(n,1)) - 1;
mu_star = ex01_bisection(f, min(x)-1, max(x)+1, 10^-10);
proj = max(x-mu_star*ones(n,1), zeros(n,1));

%% Computing the projection via cvx
cvx_begin quiet
    variable proj_cvx(n)
    minimize( 0.5*(proj_cvx-x)'*(proj_cvx-x) );
    subject to
        ones(n,1)'*proj_cvx==1;
        proj_cvx>=0;
cvx_end

%% Plots - Messages
i=0:0.1:1; %Simulate the unit simplex
j=1-i;
figure; hold on;
plot(i,j,'b-', 'DisplayName', 'Unit Simplex');
plot(x(1), x(2), 'ro', 'DisplayName', 'Point x');
plot([x(1),proj(1)], [x(2),proj(2)], 'k-', 'DisplayName', 'Distance of x from the unit simplex');
plot(proj(1), proj(2), 'ks', 'DisplayName', 'Projection of x onto the unit simplex');
legend; axis equal; grid on; xlabel('x_1'); ylabel('x_2');

fprintf('The projection onto the unit simplex computed by the corrolary formula is: [%d, %d]\n', proj);
fprintf('The projection onto the unit simplex computed by CVX is: [%d, %d]\n', proj_cvx);