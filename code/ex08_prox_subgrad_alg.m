function [x, F_out] = ex08_prox_subgrad_alg(x0, F, subgrad_f, stepsize, lambda, F_opt, c)
%   EX08_PROX_SUBGRAD_ALG Implements the proximal subgradient descent algorithm.
%
%   [X, F_OUT] = EX08_PROX_SUBGRAD_ALG(X0, F, SUBGRAD_F, STEPSIZE, LAMBDA, F_OPT, C)
%   Executes the proximal subgradient descent algorithm for function F. 
%   SUBGRAD_F is a function which belongs to the subdifferential set of F,
%   and X0 the starting point of the algorithm. Moreover, STEPSIZE function
%   determines the step size of the algorithm, and C*F_OPT is used to 
%   determine when the solution is accurate enough. Finally, LAMBDA parameter
%   is required by the algorithm. The return variables contain the points
%   X and the function values F_OUT of the algorithmic steps, from the 
%   starting point until the end.

    x = x0; %Initial point
    F_out = F(x0); %Value of f at x0
    
    k=1;
    while F_out(k) > c*F_opt  %Iterations until convergence
        x = [x, A6_soft_thresh(x(:,k)- stepsize(x(:,k),k-1)*subgrad_f(x(:,k)), lambda*stepsize(x(:,k),k-1))];
        F_out = [F_out; F(x(:,k+1))];
        k = k+1;
    end
end