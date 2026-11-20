function W = noise_directions_ou(~,Vk)
%NOISE_DIRECTIONS_OU Noise-generated directions for the OU problem.
%   W = NOISE_DIRECTIONS_OU(U,Vk) returns the additional directions
%   associated with the diffusion terms.
%
%   Inputs
%   ------
%   U  : physical basis (unused for additive noise).
%   Vk : d-by-m matrix of forced eigenvectors.
%
%   Output
%   ------
%   W : d-by-m noise-direction matrix.

    W = Vk;

end
