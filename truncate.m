function [U,Y] = truncate(X0,r)
%TRUNCATE Compress an ensemble to a prescribed rank.
%   [U,Y] = TRUNCATE(X0,r) uses the leading r left singular vectors
%   of X0, with Y = X0'*U, so that U*Y' is its rank-r approximation.
%
%   Inputs
%   ------
%   X0 : d-by-N ensemble.
%   r  : prescribed rank, 1 <= r <= min(d,N).
%
%   Outputs
%   -------
%   U : d-by-r orthonormal basis.
%   Y : N-by-r coefficient matrix.

    [U,~,~] = svd(X0,'econ');
    U = U(:,1:r);
    Y = X0'*U;

end
