function x_out = ex06_soft_thresh(x, lambda)
%   EX06_SOFT_THRESH Implements the soft thresholding operator.
%
%   X_OUT = EX06_SOFT_THRESH(X, LAMBDA) Implements the soft thresholding
%   operation, as described in Beck's book, for vector X with parameter
%   LAMBDA.

    if lambda<=0
        error('lambda should be positive');
    end
    x_out = zeros(size(x));
    x_out(abs(x)>=lambda) = sign(x(abs(x)>=lambda)).*(abs(x(abs(x)>=lambda))-lambda);
end
