function [x, f_out] = ex12_DPG(A, f, sigma, argmax_sol_func, prox_func, f_opt, epsilon)
%   EX12_DPG Implements the DPG algorithm.
%
%   [X, F_OUT] = EX12_DPG(A, F, SIGMA, ARGMAX_SOL_FUNC, PROX_FUNC, F_OPT, EPSILON)
%   Executes the DPG algorithm for function F. Functions ARGMAX_SOL_FUNC,
%   PROX_FUNC are required by the algorithm and the respective exercise, as
%   well as parameters SIGMA and matrix A. F_OPT is the optimal point of F,
%   which in combination with EPSILON form the termination criterion of the
%   algorithm. The return variables contain the points X and Z and the 
%   function values F_OUT of the algorithmic steps, from the starting point
%   until the end.

    m = size(A,1);
    y = rand(m,1);
    L = norm(A)^2/sigma;
    f_out = f_opt + 10*epsilon;
    
    k=1;
    while (abs(f_out(k) - f_opt) > epsilon) %Iterations until convergence
        x(:,k) = argmax_sol_func(y);
        y = y - A*x(:,k)/L + prox_func(A*x(:,k)-L*y)/L;
        f_out(k+1) = f(x(:,k));
        k = k+1;
    end
    
    f_out = f_out(2:end);
end