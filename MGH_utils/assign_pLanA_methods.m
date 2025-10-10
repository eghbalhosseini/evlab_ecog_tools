function pLanA_vals = assign_pLanA_methods(results, V, method, varargin)
% ASSIGN_PLANA_METHODS  Assign pLanA probabilities to depth electrodes.
% 
%   pLanA_vals = assign_pLanA_methods(results, V, method, ...)
% 
% Inputs:
%   - results: table with columns R, A, S, and optionally
%       trajectoryStart, trajectoryEnd.
%   - V:        SPM volume struct for lanANii (from spm_vol).
%   - method:   one of:
%       'trilinear', 'surface', 'radial', 'distweight', 'gaussian'
%   - varargin: additional arguments as needed:
%       * 'surfaceMesh', surf (struct with vertices Nx3) for 'surface'
%       * 'nSteps' (int) for 'radial'
%       * 'sigma' (1×3) smoothing FWHM for 'gaussian'
% 
% Output:
%   pLanA_vals: vector of length height(results)

% Precompute inverse transform
Tinv = inv(V.mat);
N = height(results);
pLanA_vals = zeros(N,1);

switch lower(method)
  case 'trilinear'
    for i=1:N
      xyz = [results.R(i); results.A(i); results.S(i); 1];
      vox = Tinv * xyz;
      pLanA_vals(i) = spm_sample_vol(V, vox(1), vox(2), vox(3), 1);
    end

  case 'surface'
    if isempty(varargin)
      error('Surface method requires a surface mesh as input');
    end
    surf = varargin{1};
    MdlKDT = KDTreeSearcher(surf.vertices);
    for i=1:N
      pt = [results.R(i), results.A(i), results.S(i)];
      idx = knnsearch(MdlKDT, pt);
      vtx = surf.vertices(idx,:);
      vox = Tinv * [vtx 1]';
      pLanA_vals(i) = spm_sample_vol(V, vox(1), vox(2), vox(3), 1);
    end

  case 'radial'
    % Check if trajectory columns exist
    if ~ismember('trajectoryStart', results.Properties.VariableNames) || ...
       ~ismember('trajectoryEnd', results.Properties.VariableNames)
      warning('trajectoryStart/trajectoryEnd not found. Using vertical sampling instead.');
      % Default to vertical sampling (5mm above and below)
      nSteps = varargin{1};
      stepSize = 1; % mm
      for i=1:N
        center = [results.R(i), results.A(i), results.S(i)];
        samples = zeros(nSteps,1);
        for k=1:nSteps
          % Sample along superior-inferior axis
          pt = center + [0, 0, (k-ceil(nSteps/2))*stepSize];
          vox = Tinv * [pt 1]';
          samples(k) = spm_sample_vol(V, vox(1), vox(2), vox(3), 1);
        end
        pLanA_vals(i) = max(samples);
      end
    else
      nSteps = varargin{1};
      for i=1:N
        origin = results.trajectoryStart(i,:);
        dirvec = results.trajectoryEnd(i,:) - origin;
        dirvec = dirvec / norm(dirvec);
        samples = zeros(nSteps,1);
        for k=1:nSteps
          pt = origin + (k-1)*1*dirvec;
          vox = Tinv * [pt 1]';
          samples(k) = spm_sample_vol(V, vox(1), vox(2), vox(3), 1);
        end
        pLanA_vals(i) = max(samples);
      end
    end

  case 'distweight'
    Y = spm_read_vols(V);
    for i=1:N
      xyz = [results.R(i); results.A(i); results.S(i); 1];
      vox = Tinv * xyz;
      i0 = floor(vox(1)); j0 = floor(vox(2)); k0 = floor(vox(3));
      offsets = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
      vals = zeros(8,1);
      dists = zeros(8,1);
      for c=1:8
        idxs = offsets(c,:) + [i0,j0,k0];
        ix=idxs(1); jy=idxs(2); kz=idxs(3);
        if ix>=1&&ix<=size(Y,1)&&jy>=1&&jy<=size(Y,2)&&kz>=1&&kz<=size(Y,3)
          vals(c)=Y(ix,jy,kz);
          dists(c)=norm([ix;jy;kz]-vox(1:3));
        else
          vals(c)=0; dists(c)=Inf;
        end
      end
      w = 1./dists; w(isinf(w))=0;
      if sum(w) == 0
        pLanA_vals(i) = 0;
      else
        pLanA_vals(i) = sum(w.*vals)/sum(w);
      end
    end

  case 'gaussian'
     % Gaussian smoothing using spm_smooth in memory
    if isempty(varargin)
        sigma = [2 2 2];  % default FWHM in mm
    else
        sigma = varargin{1};
    end
    % Read original volume
    Y = spm_read_vols(V);
    % Preallocate destination array
    Ysm = zeros(size(Y));
    % Perform smoothing in memory
    spm_smooth(Y, Ysm, sigma, 0);
    % Sample from smoothed volume
    for i = 1:N
        xyz = [results.R(i); results.A(i); results.S(i); 1];
        vox = Tinv * xyz;
        % trilinear interpolation from Ysm
        pLanA_vals(i) = interpn(Ysm, vox(1), vox(2), vox(3), 'linear', 0);
    end

  otherwise
    error('Unknown method "%s".', method);
end
end
