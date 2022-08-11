function [x, z, F_out] = ex11_ADMM(x0, z0, upd_x, upd_z, upd_y, F, F_opt, epsilon)
%   EX11_ADMM Implements the ADMM algorithm.
%
%   [X, Z, F_OUT] = EX12_ADMM(X0, Z0, UPD_X, UPD_Z, UPD_Y, F, F_OPT, EPSILON)
%   Executes the ADMM algorithm for function F, starting from points X0, Z0.
%   Functions UPD_X, UPD_Z, UPD_Y define the iterative updates of the
%   respective variables. F_OPT is the optimal point of F, which in 
%   combination with EPSILON form the termination criterion of the algorithm.
%   The return variables contain the points X and Z and the function values
%   F_OUT of the algorithmic steps, from the starting point until the end.

    y = 10*rand(length(x0),1);
    x = x0; %Initial point
    z = z0;
    F_out = F(x0, z0); %Value of F at (x0, z0)
    
    k=1;
    while abs(F_out(k)-F_opt)>epsilon  %Iterations until convergence
        x(:,k+1) = upd_x(x(:,k), z(:,k), y);
        z(:,k+1) = upd_z(x(:,k+1), z(:,k), y);
        y = upd_y(x(:,k+1), z(:,k+1), y);
        F_out(k+1) = F(x(:,k+1), z(:,k+1));
        k = k+1;
    end
end