function [x, fx] = ex07_subgrad_alg(x0, f, subgrad, stepsize, iters)
%   EX07_SUBGRAD_ALG Implements the projected subgradient descent algorithm.
%
%   [X, FX] = EX07_SUBGRAD_ALG(X0, F, SUBGRAD, STEPSIZE, ITERS)  Executes
%   the projected subgradient descent algorithm for function F. SUBGRAD is
%   a function which belongs to the subdifferential set of F, and X0 the
%   starting point of the algorithm. Moreover, STEPSIZE function determines
%   the step size of the algorithm, and ITERS determines the iterations of
%   the algorithm. The return variables contain the points X and the 
%   function values FX of the algorithmic steps, from the starting point
%   until the end.

    x = x0; %Initial point
    fx = f(x0); %Value of f at x0
    
    k=1;
    for k=1:iters %fixed number of iterations
        x = [x, x(:,k)- stepsize(x(:,k), k-1)*subgrad(x(:,k))];
        fx = [fx; f(x(:,k+1))];
    end
end