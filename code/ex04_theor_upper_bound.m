function f_best_theor = ex04_theor_upper_bound(f_opt, L_f, x0, x_opt, iters)
%   EX04_THEOR_UPPER_BOUND Provides a theoretical upper bound for the
%   Polyak step's convergence rate.
%
%   F_BEST_TH = EX04_THEOR_UPPER_BOUND(F_OPT, LF, X0, X_OPT, ITERS)
%   Provides a theoretical upper bound for the convergence rate of the
%   subgradient descent algorithm with a Polyak stepsize. F_OPT is the
%   optimal value, LF a constant described in Beck's book, X0 the starting
%   point of the algorithm, X_OPT the optimal point and ITERS the number
%   of the algorithm's iterations. Returns a vector with the function's
%   values at every iteration, according to the theoretical convergence
%   rate.

    f_best_theor = zeros(iters, 1);
    
    for k=0:iters-1
        f_best_theor(k+1) = f_opt + L_f*norm(x0-x_opt)/sqrt(k+1);
    end
end