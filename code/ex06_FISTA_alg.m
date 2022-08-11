function [x, f1_out] = ex06_FISTA_alg(x0, f1, grad_f, Lf, lambda, f1_opt, c)
%   EX06_FISTA_ALG Implements the FISTA algorithm.
%
%   [X, F1_OUT] = EX06_FISTA_ALG(X0, F1, GRAD_F, LF, LAMBDA, F1_OPT, C)
%   Executes the FISTA algorithm for function F1. GRAD_F is the gradient
%   of another function F, as described in Beck's book and X0 the starting
%   point of the algorithm. Moreover, LF is a parameter described in Beck's
%   book and depending on function F. Finally, C*F_OPT is used to determine when
%   the solution is accurate enough. The return variables contain the 
%   points X and the function values F1_OUT of the algorithmic steps, from 
%   the starting point until the end.

    x = x0; %Initial point
    f1_out = f1(x0, lambda); %Value of f at x0
    y = x0;
    t = 1; %t of (k) iteration

    k=1;
    while f1(x(:,k), lambda) > c*f1_opt %Iterations until convergence
        x = [x, A6_soft_thresh(y - grad_f(y)/Lf, lambda/Lf)];
        t_ = (1+sqrt(1+4*t^2))/2; %t of (k+1) iteration
        y = x(:,k+1) + ((t-1)/t_)*(x(:,k+1)-x(:,k));
        f1_out = [f1_out; f1(x(:,k+1), lambda)];
        k = k+1;
        t = t_;
    end
end