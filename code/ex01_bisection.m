function c = ex01_bisection(f, a, b, tolerance)
%   EX01_BISECTION  Implements the bisection algorithm.
%
%   C = EX01_BISECTION(F, A, B, TOLERANCE) Executes the bisection algorithm
%   for function F between A and B. TOLERANCE defines the precision of the
%   solution.
%
%   For more information, see <a href="bisection: 
%   web('https://en.wikipedia.org/wiki/Bisection_method')">the description of the algorithm</a>.
    
    while true
        c = (a + b)/2; %new midpoint
        if f(c)==0 || ((b-a)/2)<tolerance %solution found
            break;
        end
        if sign(f(c)) == sign(f(a))
            a = c; %new interval
        else
            b = c; %new interval
        end
    end
end