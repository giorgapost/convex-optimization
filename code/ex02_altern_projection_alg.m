function [x_k, fv] = ex02_altern_projection_alg(x0, A, b, f_v, iters)
%   EX02_ALTERN_PROJECTION_ALG  Implements the alternating projection algorithm.
%
%   [X_K, FV] = EX02_ALTERN_PROJECTION_ALG(X0, A, B, F_V, ITERS) 
%   Runs the alternating projection algorithm for function F_V and the sets
%   S1, S2 as described in the respective exercise. The starting point is X0.
%   Also, matrices A, B are parameters that define set S1. Finally, ITERS
%   determines the iterations of the algorithm, until it terminaters. The
%   return variables contain the points X_K and the function values FV of the
%   algorithmic steps, from the starting point until the end.

    n = length(x0);
    x_k = zeros(n, iters+1);
    fv = zeros(1, iters+1);
    x_k(:,1) = x0; %Initial point
    fv(1) = f_v(x0);
    
    for i=1:iters
        x_k(:,i+1) = max(x_k(:,i)-A'*((A*(A'))\(A*x_k(:,i)-b)), zeros(n,1));
        fv(i+1) = f_v(x_k(:,i+1));
    end
end