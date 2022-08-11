function [x, F_out] = ex09_SFISTA_alg(x0, F, grad_F, proj_on_set, F_opt, epsilon)
%   EX09_SFISTA_ALG Implements the S-FISTA algorithm.
%
%   [X, F_OUT] = EX09_SFISTA_ALG(X0, F, GRAD_F, PROJ_ON_SET, F_OPT, EPSILON)
%   Implements the S-FISTA algorithm for function F, starting from point
%   X0. GRAD_F is the gradient function of F. PROJ_ON_SET function projects
%   a vector on a given set, as required by the respective exercise and the
%   algorithm. F_OPT is the optimal point of F, which in combination with
%   EPSILON form the termination criterion of the algorithm. The return
%   variables contain the points X and the function values F_OUT of the
%   algorithmic steps, from the starting point until the end.

    x = x0; %Initial point
    F_out = F(x0); %Value of f at x0
    
    n = length(x);
    alpha = 1;  mu = epsilon/n;
    y = x0;    t = 1;    L =  alpha/mu;
    
    k=1;
    while F(x(:,k))-F_opt > epsilon  %Iterations until convergence
        x = [x, proj_on_set(y - grad_F(y,mu)/L)];
        t_ = (1+sqrt(1+4*t^2))/2; %t of (k+1) iteration
        y = x(:,k+1) + ((t-1)/t_)*(x(:,k+1)-x(:,k));
        F_out = [F_out; F(x(:,k+1))];
        k = k+1;
        t = t_;
    end
end