function ex05_plotting_function(f_opt, f_pol, f_dyn, f_best_th, m, n, GAUS, UNIF, MAT_TYPES)
%   EX05_PLOTTING_FUNCTION Constructs the necessary plots for ex. 5.

    for type = MAT_TYPES
        f_pol_best{type} = zeros(size(f_pol{type}));
        for k=1:length(f_pol{type})
            f_pol_best{type}(k) = min(f_pol{type}(1:k));
        end
        
        f_dyn_best{type} = zeros(size(f_dyn{type}));
        for k=1:length(f_dyn{type})
            f_dyn_best{type}(k) = min(f_dyn{type}(1:k));
        end
    end
    
    for type = MAT_TYPES
        %% Plots A5.b
        figure;
        semilogy(0:length(f_pol{type})-1, f_pol_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k_{best}-f_{opt}');
        grid on; xlim([0,length(f_pol{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question B, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question B, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A5.c
        figure;
        semilogy(0:length(f_dyn{type})-1, f_dyn_best{type}-f_opt(type), 'b-', 'DisplayName', 'f^k_{best}-f_{opt}');
        grid on; xlim([0,length(f_dyn{type})-1]); legend();
        xlabel('Iterations (k)'); ylabel('f(x_k)');
        if type==GAUS
            title(sprintf('Question C, Gaussian matrices, m=%d, n=%d',m,n));
        elseif type==UNIF
            title(sprintf('Question C, Uniform matrices, m=%d, n=%d',m,n));
        end
        
        %% Plots A5.d, A5.e
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
    end
end