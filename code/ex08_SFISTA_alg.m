function [x, F_out] = ex08_SFISTA_alg(x0, F, grad_f_smooth, Lf, alpha, beta, lambda, F_opt, epsilon)
%   EX08_SFISTA_ALG Implements the S-FISTA algorithm.
%
%   [X, F_OUT] = EX08_SFISTA_ALG(X0, F, GRAD_F_SMOOTH, LF, ALPHA, BETA, LAMBDA, F_OPT, EPSILON)
%   Implements the S-FISTA algorithm for function F, starting from point
%   X0. GRAD_F_SMOOTH is the gradient function of the smoothened F. LF, 
%   ALPHA, BETA, LAMBDA are defined in Beck's book and are related to the
%   S-FISTA algorithm and the smoothing process. F_OPT is the optimal point
%   of F, which in combination with EPSILON form the termination criterion
%   of the algorithm. The return variables contain the points X and the
%   function values F_OUT of the algorithmic steps, from the starting point
%   until the end.

    x = x0; %Initial point
    F_out = F(x0); %Value of f at x0
    
    mu = sqrt(alpha/beta)*epsilon / (sqrt(alpha*beta) + sqrt(alpha*beta + Lf*epsilon));
    L =  Lf + alpha/mu;

    y = x0;
    t = 1; %t of (k) iteration
        
    k=1;
    while F(x(:,k))-F_opt > epsilon  %Iterations until convergence
        x = [x, A6_soft_thresh(y - grad_f_smooth(y, mu)/L, lambda/L)];
        t_ = (1+sqrt(1+4*t^2))/2; %t of (k+1) iteration
        y = x(:,k+1) + ((t-1)/t_)*(x(:,k+1)-x(:,k));
        F_out = [F_out; F(x(:,k+1))];
        k = k+1;
        t = t_;
    end
end