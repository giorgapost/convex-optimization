function [x, F_out] = ex13_FDPG(n, F, L, D, argmax_sol_func, prox_func, F_opt, epsilon)
%   EX13_FDPG Implements the DPG algorithm.
%
%   [X, F_OUT] = EX13_FDPG(N, F, L, D, ARGMAX_SOL_FUNC, PROX_FUNC, F_OPT, EPSILON)
%   Executes the FDPG algorithm for function F. Functions ARGMAX_SOL_FUNC,
%   PROX_FUNC, D are required by the algorithm and the respective exercise,
%   as well as parameters N,L. F_OPT is the optimal point of F, which in 
%   combination with EPSILON form the termination criterion of the
%   algorithm. The return variables contain the points X and Z and the 
%   function values F_OUT of the algorithmic steps, from the starting point
%   until the end.

    y=rand(n-1,1); w=y; t=1;
    F_out = F_opt + 10*epsilon;
    
    k=1;
    while (abs(F_out(k) - F_opt) > epsilon) %Iterations until convergence
        u = argmax_sol_func(w);
        y_ = w - D(u)/L + prox_func(D(u)-L*w)/L;
        t_ = (1+sqrt(1+4*t^2))/2;
        w = y_ + ((t-1)/t_)*(y_-y);
        x(:,k) = argmax_sol_func(y_);
        
        F_out(k+1) = F(x(:,k));
        k = k+1;
        y = y_;
        t = t_;
    end
    
    F_out = F_out(2:end);
end