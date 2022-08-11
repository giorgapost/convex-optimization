function c = ex10_bisection(f, x, a, b, tolerance)
%   EX10_BISECTION  Implements the bisection algorithm.
%
%   C = EX10_BISECTION(F, X, A, B, TOLERANCE) Executes the bisection algorithm
%   for function F(X,Y) between A and B, with respect to Y.
%   TOLERANCE defines the precision of the solution.
%
%   For more information, see <a href="bisection: 
%   web('https://en.wikipedia.org/wiki/Bisection_method')">the description of the algorithm</a>.

    while true
        c = (a + b)/2; %new midpoint
        if f(x,c)==0 || ((b-a)/2)<tolerance %solution found
            break;
        end
        if sign(f(x,c)) == sign(f(x,a))
            a = c; %new interval
        else
            b = c; %new interval
        end
    end
end