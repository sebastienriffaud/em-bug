function error = strong_error(Xexact,X)
%STRONG_ERROR Compute the maximum-in-time RMS strong ensemble error.
%   error = STRONG_ERROR(Xexact,X) returns
%       max_k norm(Xexact(:,:,k)-X(:,:,k),'fro')/sqrt(N).
%
%   Inputs
%   ------
%   Xexact : d-by-N-by-(Nt+1) exact trajectory.
%   X      : d-by-N-by-(Nt+1) approximate trajectory.
%
%   Output
%   ------
%   error : nonnegative scalar strong error.

    difference = reshape(Xexact-X,[],size(X,3));
    error = max(vecnorm(difference,2,1))/sqrt(size(X,2));

end
