function [x, F_out] = ex07_SFISTA_alg(x0, F, grad_f, Lf, D, D_T, normD, F_opt, epsilon)
%   EX07_SFISTA_ALG Implements the S-FISTA algorithm.
%
%   [X, F_OUT] = EX07_SFISTA_ALG(X0, F, GRAD_F, LF, LAMBDA, D, D_T, NORMD, F_OPT, EPSILON)
%   Implements the S-FISTA algorithm for function F, starting from point
%   X0. GRAD_F is the gradient function of F. LF, D, D_T, normD are defined
%   in Beck's book. F_OPT is the optimal point of F, which in combination
%   with EPSILON form the termination criterion of the algorithm. The return
%   variables contain the points X and the function values F_OUT of the
%   algorithmic steps, from the starting point until the end.

    x = x0; %Initial point
    F_out = F(x0); %Value of f at x0
    
    n = length(x);
    alpha = normD^2; beta = n/2;
    mu = sqrt(alpha/beta)*epsilon / (sqrt(alpha*beta) + sqrt(alpha*beta + Lf*epsilon));
    
    y = x0;
    t = 1; %t of (k) iteration
    L =  Lf + alpha/mu;
    
    grad_F = @(x) grad_f(x) + (1/mu)*D_T((D(x)-A6_soft_thresh(D(x), mu)));
    
    k=1;
    while F(x(:,k))-F_opt > epsilon  %Iterations until convergence
        x = [x, y - grad_F(y)/L];
        t_ = (1+sqrt(1+4*t^2))/2; %t of (k+1) iteration
        y = x(:,k+1) + ((t-1)/t_)*(x(:,k+1)-x(:,k));
        F_out = [F_out; F(x(:,k+1))];
        k = k+1;
        t = t_;
    end
end