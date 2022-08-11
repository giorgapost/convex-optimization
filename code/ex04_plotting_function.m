function ex04_plotting_function(f_opt, f_pol, f_dyn, f_best_th, f_stoch, f_incr, m, n, GAUS, UNIF, MAT_TYPES)
%   EX04_PLOTTING_FUNCTION Constructs the necessary plots for ex. 4.
    
    for type = MAT_TYPES
        f_pol_best{type} = zeros(size(f_pol{type}));
        for k=1:length(f_pol{type})
            f_pol_best{type}(k) = min(f_pol{type}(1:k));
        end
        
        f_dyn_best{type} = zeros(size(f_dyn{type}));
        for k=1:length(f_dyn{type})
            f_dyn_best{type}(k) = min(f_dyn{type}(1:k));
        end
        
        f_stoch_best{type} = zeros(size(f_stoch{type}));
        for k=1:length(f_stoch{type})
            f_stoch_best{type}(k) = min(f_stoch{type}(1:k));
        end
        
        f_incr_best{type} = zeros(size(f_incr{type}));
        for k=1:length(f_incr{type})
            f_incr_best{type}(k) = min(f_incr{type}(1:k));
        end
    end
    
    for type = MAT_TYPES
        %% Plots A4.b
        figure;
        semilogy(0:length(f_pol{type})-1, f_pol_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k_{best}-f_{opt}');
        grid on; xlim([0,length(f_pol{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question B, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question B, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A4.c
        figure;
        semilogy(0:length(f_dyn{type})-1, f_dyn_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k_{best}-f_{opt}');
        grid on; xlim([0,length(f_dyn{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question C, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question C, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A4.d, A4.e
        figure;
        semilogy(0:length(f_pol{type})-1, f_pol_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^{k(Polyak)}_{best}-f_{opt}');
        hold on;
        semilogy(0:length(f_dyn{type})-1, f_dyn_best{type}-f_opt(type), 'r-', 'DisplayName', 'f^{k(dynamic)}_{best}-f_{opt}');
        semilogy(0:length(f_best_th{type})-1, f_best_th{type}-f_opt(type), 'k-', 'DisplayName', 'Upper bound for f^{k}_{best}-f_{opt}');
        grid on; xlim([0,max(length(f_dyn{type})-1, length(f_pol{type})-1)]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Questions D,E, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Questions D,E, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A4.f
        figure;
        semilogy(0:length(f_stoch{type})-1, f_stoch{type}-f_opt(type), 'r-', 'DisplayName', 'f(x_{km})-f_{opt}');
        hold on;
        semilogy(0:length(f_stoch_best{type})-1, f_stoch_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k^m_{best}-f_{opt}');
        grid on; xlim([0,length(f_stoch{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k_m)');
        if type==GAUS
            title(sprintf('Question F, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question F, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A4.g
        figure;
        semilogy(0:length(f_incr{type})-1, f_incr{type}-f_opt(type), 'r-', 'DisplayName', 'f(x_{k})-f_{opt}');
        hold on;
        semilogy(0:length(f_incr_best{type})-1, f_incr_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k_{best}-f_{opt}');
        grid on; xlim([0,length(f_incr{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question G, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question G, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A4.h
        figure;
        semilogy(0:length(f_pol{type})-1, f_pol_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^{k(Polyak step)}_{best}-f_{opt}');
        hold on;
        semilogy(0:length(f_dyn{type})-1, f_dyn_best{type}-f_opt(type), 'm-', 'DisplayName', 'f^{k(dynamic step)}_{best}-f_{opt}');
        semilogy(0:length(f_best_th{type})-1, f_best_th{type}-f_opt(type), 'k-', 'DisplayName', 'Upper bound for f^{k}_{best}-f_{opt}');
        semilogy(0:length(f_stoch{type})-1, f_stoch_best{type}-f_opt(type), 'g-', 'DisplayName', 'f^{km(stochastic subgr.)}_{best}-f_{opt}');
        semilogy(0:length(f_incr{type})-1, f_incr_best{type}-f_opt(type), 'r-', 'DisplayName', 'f^{k(incremental subgr.)}_{best}-f_{opt}');
        grid on; xlim([0,max(length(f_dyn{type})-1, length(f_pol{type})-1)]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question H, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question H, Uniform matrices, m=%d, n=%d',m,n));
        end
    end
end