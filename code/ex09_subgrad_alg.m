function [x, f_out] = ex09_subgrad_alg(x0, f, subgrad, stepsize, proj_on_set, iters)
%   EX09_SUBGRAD_ALG Implements the projected subgradient descent algorithm.
%
%   [X, F_OUT] = EX09_SUBGRAD_ALG(X0, F, SUBGRAD, STEPSIZE, PROJ_ON_SET, ITERS)
%   Executes the projected subgradient descent algorithm for function F. 
%   SUBGRAD is a function which belongs to the subdifferential set of F,
%   and X0 the starting point of the algorithm. Moreover, STEPSIZE function
%   determines the step size of the algorithm, and ITERS determines the 
%   iterations of the algorithm. Also, PROJ_ON_SET function projects
%   a vector on a specific set. The return variables contain the points X
%   and the function values FX of the algorithmic steps, from the starting
%   point until the end.

    x = x0; %Initial point
    f_out = f(x0); %Value of f at x0
    
    k=1;
    for k=1:iters %fixed number of iterations
        x = [x, proj_on_set(x(:,k)- stepsize(x(:,k), k-1)*subgrad(x(:,k)))];
        f_out = [f_out; f(x(:,k+1))];
    end
end