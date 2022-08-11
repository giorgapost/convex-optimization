function [x, f1_out] = ex06_ISTA_alg(x0, f1, grad_f, Lf, lambda, f1_opt, c)
%   EX06_ISTA_ALG Implements the ISTA algorithm.
%
%   [X, F1_OUT] = EX06_ISTA_ALG(X0, F1, GRAD_F, LF, LAMBDA, F1_OPT, C)
%   Executes the ISTA algorithm for function F1. GRAD_F is the gradient
%   of another function F, as described in Beck's book and X0 the starting
%   point of the algorithm. Moreover, LF is a parameter described in Beck's
%   book and depending on function F. Finally, F1_OPT is the optimal point
%   of F1 which, in combination with C, determines the termination criterion
%   of the algorithm. The return variables contain the points X and the
%   function values F1_OUT of the algorithmic steps, from the starting point
%   until the end.

    x = x0; %Initial point
    f1_out = f1(x0, lambda); %Value of f at x0

    k=1;
    while f1(x(:,k), lambda) > c*f1_opt %Iterations until convergence
        x = [x, A6_soft_thresh(x(:,k) - grad_f(x(:,k))/Lf, lambda/Lf)];
        f1_out = [f1_out; f1(x(:,k+1), lambda)];
        k = k+1;
        k
    end
end