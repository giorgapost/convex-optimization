function [x, fx] = ex04_subgrad_stochastic_alg(x0, f, subgrad_fi, stepsize, epochs, m)
%   EX04_SUBGRAD_STOCHASTIC_ALG Implements the stochastic subgradient descent algorithm.
%
%   [X, FX] = EX04_SUBGRAD_STOCHASTIC_ALG(X0, F, SUBGRAD_FI, STEPSIZE, EPOCHS, M)
%   Executes the stochastic subgradient descent algorithm for function F. 
%   SUBGRAD_FI is a function which belongs to the subdifferential set of Fi,
%   i=1,...,M and X0 the starting point of the algorithm. Moreover, STEPSIZE
%   function determines the step size of the algorithm, and EPOCHS is the
%   number of iterations for the algorithm. The return variables contain the
%   points X and the function values FX of the algorithmic steps, from the
%   starting point until the end.

    x = x0; %Initial point
    fx = f(x0); %Value of f at x0
    
    x_iter = x0;
    for e=1:epochs
        for mm = 1:m
            k = (e-1)*m+mm-1;
            i = randi([1,m]);
            x_iter = x_iter - stepsize(k)*subgrad_fi(i, x_iter);
        end
        x = [x, x_iter];
        fx = [fx; f(x_iter)];
    end
end