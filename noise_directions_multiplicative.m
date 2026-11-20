function W = noise_directions_multiplicative(U,V,mu)
%NOISE_DIRECTIONS_MULTIPLICATIVE Noise-generated directions for the multiplicative problem.
%   W = NOISE_DIRECTIONS_MULTIPLICATIVE(U,V,mu) returns the additional
%   directions associated with the diffusion terms.
%
%   Inputs
%   ------
%   U  : d-by-r physical basis.
%   V  : d-by-d orthonormal modal basis.
%   mu : d-by-m diffusion coefficients.
%
%   Output
%   ------
%   W : d-by-(r*m) noise-direction matrix.

    Z = V'*U;
    W = zeros(size(V,1),size(U,2),size(mu,2));
    for j = 1:size(mu,2)
        W(:,:,j) = V*(mu(:,j).*Z);
    end
    W = reshape(W,size(V,1),[]);

end
