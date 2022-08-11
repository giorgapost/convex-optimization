% Author: Georgios Apostolakis
close all;
clear;
clc;

%% Definition of constants & functions, Setup of the problem
GAUS = 1; %Gaussian distribution
UNIF = 2; %Uniform distribution
MAT_TYPES = [GAUS, UNIF];
DIMENSIONS = [ 6, 5;  %(m, n) dimensions of the matrices
              25, 5];
rng(0);
for dim_ptr = 1:size(DIMENSIONS,1)
    m = DIMENSIONS(dim_ptr,1);
    n = DIMENSIONS(dim_ptr,2);

%% Initialization
    c = (dim_ptr==2)*1.01 + (dim_ptr~=2)*1.1; %Termination criterion    
    x0 = rand(n,1);  %Initial point
    x_opt = zeros(n, length(MAT_TYPES));
    f_opt = zeros(length(MAT_TYPES), 1);
    A = zeros(m,n,length(MAT_TYPES));
    b = zeros(m,length(MAT_TYPES));
    for type = MAT_TYPES
        if type==GAUS
            A(:,:,type) = 10*randn(m,n);  %Gaussian generated
            b(:,type) = 10*randn(m,1);  %Gaussian generated
        elseif type==UNIF
            A(:,:,type) = 10*rand(m,n);  %Uniformly generated
            b(:,type) = 10*rand(m,1);  %Uniformly generated
        end
    end
    f = @(x, type) norm(A(:,:,type)*x-b(:,type), 1);  %Cost function
    sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book
    subgrad_f = @(x, type) A(:,:,type)'*sgn(A(:,:,type)*x-b(:,type));  %subgradient from the subdifferential

    for type = MAT_TYPES
        cvx_begin quiet
            variable opt_pnt(n)
            minimize f(opt_pnt, type)
        cvx_end
        x_opt(:,type) = opt_pnt;
        f_opt(type) = cvx_optval;  %f_opt for the Polyak step and the gaussian matrices
    end

%% Solving with the projected subgradient descent algorithm
    for type = MAT_TYPES
        ff = @(x) f(x,type);
        subgrad_ff = @(x) subgrad_f(x,type);
        
        %4b. Subgradient algorithm with Polyak step
        polyak_step = @(x, k) (norm(subgrad_ff(x))>0).*(ff(x)-f_opt(type)) / (norm( subgrad_ff(x) )^2) + (norm(subgrad_ff(x))==0)*1;
        [x_pol{type}, f_pol{type}] = ex04_subgrad_alg(x0, ff, subgrad_ff, polyak_step, f_opt(type), c);
    
        %4c. Subgradient algorithm with Dynamic step
        L_f = norm(A(:,:,type))*sqrt(m);
        dynamic_step = @(x, k) (norm(subgrad_ff(x))>0).*1/(norm(subgrad_ff(x))*sqrt(k+1)) + (norm(subgrad_ff(x))==0)*(1/L_f);
        [x_dyn{type}, f_dyn{type}] = ex04_subgrad_alg(x0, ff, subgrad_ff, dynamic_step, f_opt(type), c);
        
        %4e. Theoretical upper bound
        f_best_th{type} = ex04_theor_upper_bound(f_opt(type), L_f, x0, x_opt(:,type), max(length(f_dyn{type}),length(f_pol{type})));
        
        %4f. Stochastic subgradient algorithm
        theta = norm(x0-x_opt(:,type))^2;
        L_f = sqrt(m)*norm(A(:,:,type));
        stepsize = @(k) sqrt(2*theta)*m /(L_f*sqrt(k+1));
        subgrad_fi = @(i,x) (A(i,:,type)*x-b(i,type)>=0)*(A(i,:,type)') + (A(i,:,type)*x-b(i,type)<0)*(-A(i,:,type)');
        [x_stoch{type}, f_stoch{type}] = ex04_subgrad_stochastic_alg(x0, ff, subgrad_fi, stepsize, length(f_dyn{type})-1, m);
        
        %4g. Incremental subgradient algorithm
        [x_incr{type}, f_incr{type}] = ex04_subgrad_stochastic_alg(x0, ff, subgrad_fi, stepsize, length(f_dyn{type})-1, m);
    end
    
%% Plots
    ex04_plotting_function(f_opt, f_pol, f_dyn, f_best_th, f_stoch, f_incr, m, n, GAUS, UNIF, MAT_TYPES)
end