function subgrad_f = ex05_subgrad_f(x,A,b)
%   EX05_SUBGRAD_F Provides the subgradient of the function in ex. 5.
%
%   SUBGRAD_F = EX05_SUBGRAD_F(X,A,B) Returns the subgradient of the
%   function f given in ex. 5, for the provided values A,B,X of matrix A
%   and vectors b,x respectively.

    sgn = @(x) double(x>=0) - double(x<0); %as defined in Beck's book

    subgrad_f = zeros(size(x));
    i = 1;
    if sum(abs(A*x-b))~=0
        while abs(abs(A(i,:)*x-b(i)) - norm(A*x-b, Inf))>10^(-10)
            i = i+1;
        end
        subgrad_f = sgn(A(i,:)*x-b(i))*(A(i,:)');
    end
end