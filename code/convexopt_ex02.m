% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Initialization
rng(0);
m = 4;
n = 8;  %dimension of the problem
A = 10*randn(m,n);
b = 10*rand(m,1);
x0 = zeros(n,1); %>=0, belongs in S2
iters = 2500; %iterations
f_v =  @(x) max([abs(A*x-b)./vecnorm(A')'; norm(x - max(x,zeros(n,1)))]); %constraint violation function

%% Checking for a non-empty intersection of the sets
while true
   cvx_begin quiet
      variable x_feas(n)
      minimize 1
      subject to
         A*x_feas==b;
         x_feas>=0;
   cvx_end
   if ~isnan(sum(x_feas))
      fprintf("The intersection of S1, S2 is not an empty set.\n");
      break
   end
   A = rand(m,n); %If we reach here, the intersection was found to be empty
   b = rand(m,1);
end

%% Solving with the alternating projection algorithm
[~, fv_al] = ex02_altern_projection_alg(x0, A, b, f_v, iters);

%% Solving with the greedy projection algorithm
[~, fv_gr] = ex02_greedy_projection_alg(x0, A, b, f_v, iters);

%% Plots
figure;
semilogy(0:iters, fv_al, 'r-', 'DisplayName', 'Alternating projection algorithm');
hold on;
semilogy(0:iters, fv_gr, 'b-', 'DisplayName', 'Greedy projection algorithm');
legend; grid on;
xlabel('Iterations (k)'); ylabel('f_v(x_k)');