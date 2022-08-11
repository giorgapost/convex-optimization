function [x, f_out] = ex10_mirror_descent(x0, f, subgrad, stepsize, f_opt, c)
%   EX10_MIRROR_DESCENT Implements the mirror descent algorithm.
%
%   [X, F_OUT] = EX10_MIRROR_DESCENT(X0, F, SUBGRAD, STEPSIZE, F_OPT, C)
%   Executes the mirror descent algorithm for function F. SUBGRAD is
%   a function which belongs to the subdifferential set of F, and X0 the
%   starting point of the algorithm. Moreover, STEPSIZE function determines
%   the step size of the algorithm, and C*F_OPT is used to determine when
%   the solution is accurate enough. The return variables contain the points
%   X and the function values F_OUT of the algorithmic steps, from the 
%   starting point until the end.

    x = x0; %Initial point
    f_out = f(x0); %Value of f at x0
    
    k=1;
    while f_out(k) > c*f_opt  %Iterations until convergence
        sum0 = sum(x(:,k).*exp(-stepsize(x(:,k),k-1)*subgrad(x(:,k))));
        x_next = x(:,k).*exp(-stepsize(x(:,k),k-1)*subgrad(x(:,k)))./sum0;
        x = [x, x_next];
        f_out = [f_out; f(x(:,k+1))];
        k = k+1;
    end
end