function [x, F_out] = ex13_DPG(n, F, L, D, argmax_sol_func, prox_func, F_opt, epsilon)
%   EX13_DPG Implements the DPG algorithm.
%
%   [X, F_OUT] = EX13_DPG(N, F, L, D, ARGMAX_SOL_FUNC, PROX_FUNC, F_OPT, EPSILON)
%   Executes the DPG algorithm for function F. Functions ARGMAX_SOL_FUNC,
%   PROX_FUNC, D are required by the algorithm and the respective exercise,
%   as well as parameters N,L. F_OPT is the optimal point of F, which in 
%   combination with EPSILON form the termination criterion of the
%   algorithm. The return variables contain the points X and Z and the 
%   function values F_OUT of the algorithmic steps, from the starting point
%   until the end.

    y = rand(n-1,1);
    F_out = F_opt + 10*epsilon;
    
    k=1;
    while (abs(F_out(k) - F_opt) > epsilon) %Iterations until convergence
        x(:,k) = argmax_sol_func(y);
        y = y - D(x(:,k))/L + prox_func(D(x(:,k))-L*y)/L;
        F_out(k+1) = F(x(:,k));
        k = k+1;
    end
    
    F_out = F_out(2:end);
end