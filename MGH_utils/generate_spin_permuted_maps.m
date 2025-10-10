function spin_maps = generate_spin_permuted_maps(coords, values, nPerm)
    % Generate spin-permuted maps preserving spatial autocorrelation
    %
    % Inputs:
    %   coords - Nx3 electrode coordinates (R,A,S)
    %   values - Nx1 values to permute (e.g., LanA probabilities)
    %   nPerm  - number of permutations
    %
    % Output:
    %   spin_maps - N×nPerm matrix of spin-permuted values

    N = size(coords,1);
    spin_maps = zeros(N, nPerm);

    % Project coords to unit sphere
    norms = sqrt(sum(coords.^2, 2));
    xyz_sph = coords ./ norms;

    % Build KD-tree once for nearest-neighbor lookups
    Mdl = KDTreeSearcher(xyz_sph);

    for p = 1:nPerm
        % Generate random rotation: use uniformly-sampled Euler angles
        % using the algorithm: https://github.com/frantisekvasa/rotate_parcellation
        u = rand();
        v = rand();
        w = rand();
        alpha = 2*pi*u;
        beta  = acos(2*v - 1);
        gamma = 2*pi*w;

        % Rotation matrix from Euler Z–X'–Z'' convention
        Rz1 = [ cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1 ];
        Rx  = [ 1 0 0; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta) ];
        Rz2 = [ cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1 ];
        R   = Rz2 * Rx * Rz1;

        % Rotate all spherical coordinates
        xyz_rot = (R * xyz_sph')';

        % Nearest neighbor mapping
        idx = knnsearch(Mdl, xyz_rot);

        % Assign permuted values
        spin_maps(:,p) = values(idx);
    end
end
