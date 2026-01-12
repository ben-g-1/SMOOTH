function stat = SMOOTHstat(cfg, varargin)

% SMOOTHstat - Surrogate map-based group-level analysis of human intracranial EEG data
%
% Use as:
%   stat = SMOOTHstat(cfg, source1, source2, ...)
%
% where source1..N are data structures with fields
%   source.stat         = dependent variable per channel
%   source.label        = channel labels
%   source.elec         = electrode struct with .chanpos in fsavg coords
%
% cfg fields
%   cfg.fshome           = FreeSurfer home (required)
%   cfg.parameter        = field with DV (default 'stat')
%   cfg.smooth           = 'sphere' (default) | 'pial'
%   cfg.kernel           = 'gaussian' (default) | 'boxcar'
%   cfg.kernelwidth      = FWHM (gaussian) or radius (boxcar) in mm (default 20)
%   cfg.graphsigma       = Laplacian kernel σ for spectral surrogates in mm (default 20)
%   cfg.rankrescale      = 'no' | 'quantile' | 'zscore'
%   cfg.equalize_n       = 'no' | 'yes'  (downweight high coverage in t)
%   cfg.normalize        = 'no' | 'yes'  (per-subject z before t)
%   cfg.randmethod       = 'signflip' (default) | 'bandrotate'
%   cfg.numrandomization = permutations (default 1000)
%   cfg.minnbsub         = mininum number of subjects per vertex (default 3)
%   cfg.clusteralpha     = vertex-level cluster α per tail (default 0.05)
%   cfg.keepmaps         = 'no' | 'yes'
%   cfg.keepsurrogates   = 'no' | 'yes' (memory expensive)
%
% Output (stat)
%   stat.stat, stat.prob, stat.mask, clusters (+ distributions), coverage, cfg
%
% Requirements: FreeSurfer fsaverage, FieldTrip
% Authors: Kamren Khan & Arjen Stolk, v0.1


% ---- Defaults
cfg                  = ft_checkconfig(cfg, 'required', {'fshome'});
cfg.parameter        = ft_getopt(cfg, 'parameter',         'stat');
cfg.smooth           = ft_getopt(cfg, 'smooth',          'sphere');
cfg.kernel           = ft_getopt(cfg, 'kernel',        'gaussian');
cfg.kernelwidth      = ft_getopt(cfg, 'kernelwidth',           20); % ideally informed by our fwhm analysis
cfg.graphsigma       = ft_getopt(cfg, 'graphsigma',            20);
cfg.rankrescale      = ft_getopt(cfg, 'rankrescale',         'no');
cfg.equalize_n       = ft_getopt(cfg, 'equalize_n',          'no');
cfg.normalize        = ft_getopt(cfg, 'normalize',           'no');
cfg.randmethod       = ft_getopt(cfg, 'randmethod',    'signflip');
cfg.numrandomization = ft_getopt(cfg, 'numrandomization',    1000);
cfg.minnbsub         = ft_getopt(cfg, 'minnbsub',               3); % This seems low. Maybe there should be a warning if there are fewer than 10 subjects or so? BG
cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',        0.05);
cfg.keepmaps         = ft_getopt(cfg, 'keepmaps',            'no');
cfg.keepsurrogates   = ft_getopt(cfg, 'keepsurrogates',      'no');

% ---- Meshes (fsaverage) and adjacency (for clustering)
lh      = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','lh.pial'));
rh      = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','rh.pial'));
sph_l   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','lh.sphere.reg'));
sph_r   = ft_read_headshape(fullfile(cfg.fshome,'subjects','fsaverage','surf','rh.sphere.reg'));
sph.pos = [sph_l.pos; sph_r.pos];
sph.tri = [sph_l.tri; sph_r.tri + size(sph_l.pos,1)];
i = [sph.tri(:,1); sph.tri(:,1); sph.tri(:,2); sph.tri(:,2); sph.tri(:,3); sph.tri(:,3)];
j = [sph.tri(:,2); sph.tri(:,3); sph.tri(:,1); sph.tri(:,3); sph.tri(:,1); sph.tri(:,2)];
adj = sparse(i,j,1,size(sph.pos,1),size(sph.pos,1)); adj = adj | adj'; % undirected

% -------------------------------------------------------------------------
% Step 1: Generation of Cortical Activation Maps
% -------------------------------------------------------------------------
fprintf('generating empirical maps\n')

nSubjects = numel(varargin);
nVertL  = size(sph_l.pos,1);
nVertR  = size(sph_r.pos,1);
allmaps = nan(nVertL+nVertR, nSubjects);

% convert FWHM (mm) to σ (mm) for Gaussian kernel (if used)
if strcmpi(cfg.kernel, 'gaussian')
  sigma_mm = cfg.kernelwidth/(2*sqrt(2*log(2)));   % FWHM -> σ
  % Might be helpful to add a check here to ensure that the kernel width is not too large relative to the mesh size. BG
  % Also verbosity could be added to inform the user of the sigma value being used- since gaussian is default, user might be unaware of transformation. BG
elseif strcmpi(cfg.kernel, 'boxcar')
  sigma_mm = cfg.kernelwidth;                      % boxcar radius
else
  error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
end

% one cache reused across subjects and hemispheres
e2m_cache = struct('pial', struct('L', [], 'R', []), ...  % I might refactor this into something more readable, like global_cache or similar. BG
  'sphere', struct('kdtL', [], 'kdtR', []));              % The nested structure is maybe unnecessary? Definitely unintuitive at this point, but I'll delete this comment if it makes sense later. BG

for s = 1:nSubjects
  data = varargin{s};
  dv   = data.(cfg.parameter);
  pos  = data.elec.chanpos;
  hemi = pos(:,1) < 0;  % fsavg: L<0  I would refactor so that hemi is left_hemi or hemiL following later convention for clarity. BG

  % electrode to vertex smoothing
  [mapL, cacheL, idxL, e2m_cache] = elec2map_build('L', cfg, lh, sph_l, pos, dv, sigma_mm,  hemi, e2m_cache);
  [mapR, cacheR, idxR, e2m_cache] = elec2map_build('R', cfg, rh, sph_r, pos, dv, sigma_mm, ~hemi, e2m_cache);

  allmaps(1:nVertL,               s) = mapL;
  allmaps(nVertL+1:nVertL+nVertR, s) = mapR;

  % % FIXME remove: check smoothing
  % figure;
  % hold on; ft_plot_mesh(lh, 'vertexcolor', 'curv');
  % hold on; ft_plot_mesh(lh, 'vertexcolor', mapL);
  % view([-90 0])

  % rank-rescale prep
  rescale_cache = struct();
  mL = isfinite(mapL); % refactor to avoid using same variable twice BG
  rescale_cache.L = struct('mask',mL,'tgt_sorted',sort(mapL(mL),'ascend'), ...
    'mu',mean(mapL(mL)), 'sd',max(std(mapL(mL)),eps));
  mR = isfinite(mapR);
  rescale_cache.R = struct('mask',mR,'tgt_sorted',sort(mapR(mR),'ascend'), ...
    'mu',mean(mapR(mR)), 'sd',max(std(mapR(mR)),eps));

  % keep for permutations
  data.elec.hemi     = hemi;
  data.elec.idxL     = idxL;
  data.elec.idxR     = idxR;
  data.smooth_cache  = struct('L',cacheL,'R',cacheR);
  data.rescale_cache = rescale_cache;
  varargin{s} = data;
end

% -------------------------------------------------------------------------
% Step 2: Surrogate Map Generation via Spectral Resampling
% -------------------------------------------------------------------------
fprintf('generating surrogate maps\n')

for s = 1:nSubjects
  data = varargin{s};

  % electrode sphere coords (arbitrary units) for eigenbasis
  Lh_sph = sph_l.pos(data.elec.idxL, :);
  Rh_sph = sph_r.pos(data.elec.idxR, :);

  % (a)+(b): Graph construction + Eigenbasis (per hemisphere)
  [QL, lamL] = graph_eigenbasis(Lh_sph, mm2sphere(cfg.graphsigma));
  [QR, lamR] = graph_eigenbasis(Rh_sph, mm2sphere(cfg.graphsigma));
  Q = blkdiag(QL, QR); % (nL+nR) x (nL+nR)

  % (c) Projection (modal coefficients)
  data.modes       = Q;                         % elecs x modes
  data.modalcoeffs = Q' * data.(cfg.parameter); % modal coefficients
  data.eigvals     = [lamL; lamR];              % eigenvalues
  data.B           = data.modes * spdiags(data.modalcoeffs,0,length(data.modalcoeffs),length(data.modalcoeffs));
  varargin{s} = data;
end

% -------------------------------------------------------------------------
% Steps 3: Group-Level Cluster Statistics
% -------------------------------------------------------------------------
fprintf('performing cluster-based permutation testing\n')

coverage = sum(isfinite(allmaps), 2);  % V × 1
mask     = coverage < cfg.minnbsub; % vertices below min subjects excluded
                                    % This could be made more verbose to inform the user of how many vertices were excluded based on coverage. BG
                                    % fprintf('Excluding %d vertices with coverage below %d subjects.\n', sum(mask), cfg.minnbsub);

negtailcritval = norminv(cfg.clusteralpha);
negdistribution = zeros(cfg.numrandomization,1);

postailcritval = norminv(1 - cfg.clusteralpha);
posdistribution = zeros(cfg.numrandomization,1); % Mixed order bothered me... BG

if strcmp(cfg.keepsurrogates, 'yes')
  surr_subj  = zeros(nSubjects, size(sph.pos,1), cfg.numrandomization);
  surr_group = zeros(size(sph.pos,1), cfg.numrandomization);
end

% main loop: permutations + one empirical pass
for iter = 1:cfg.numrandomization+1
  if iter <= cfg.numrandomization
    fprintf(' permutation %d/%d\n', iter, cfg.numrandomization);

    permMaps    = nan(size(allmaps)); % V × S
    % sigma5_mm   = 5/(2*sqrt(2*log(2))); % FWHM(=5mm) -> sigma(mm)
    cfg5        = cfg; 
    cfg5.kernel = 'gaussian';
    cfg5.smooth = 'sphere'; % sphere for speed (geodesic full-vertex would be too slow)

    % [~, rankSmoothL, ~, e2m_cache] = elec2map_build('L', cfg5, lh, sph_l, lh.pos, zeros(nVertL,1), sigma5_mm, true(nVertL,1), e2m_cache);
    % [~, rankSmoothR, ~, e2m_cache] = elec2map_build('R', cfg5, rh, sph_r, rh.pos, zeros(nVertR,1), sigma5_mm, true(nVertR,1), e2m_cache);

    for s = 1:nSubjects
      data = varargin{s};

      % (d) Resampling: spectral surrogate in electrode space
      switch lower(cfg.randmethod)
        case 'signflip' % random sign flips
          sflip = sign(randn(size(data.modalcoeffs)));
          dvn   = data.B * sflip;   % fast: B*sflip
        case 'bandrotate' % band-oriented orthogonal rotation
          a2    = spectral_bandrotate(data.modalcoeffs, data.eigvals);
          dvn   = data.modes * a2;  % reconstruct
        otherwise
          error('cfg.randmethod must be ''signflip'' or ''bandrotate''.');
      end

% %%%%%%%%%%
%       % RESCALE DV AT ELECTRODE LEVEL: JUST FOR TESTING
%       [~,dvidx] = sort(dvn); dvn(dvidx) = sort(data.stat);
% %%%%%%%%%%

      % (e) Postprocessing: apply same smoothing as empirical
      mapL = elec2map_apply(data.smooth_cache.L, dvn(data.elec.hemi));
      mapR = elec2map_apply(data.smooth_cache.R, dvn(~data.elec.hemi));

      % optional marginal matching
      switch lower(cfg.rankrescale)
        case 'no'
        case 'exact'
          [~, idxL] = sort(mapL); mapL(idxL) = sort(allmaps(1:163842,s));  % why is this hardcoded? Throws errors with the demo script. BG
          [~, idxR] = sort(mapR); mapR(idxR) = sort(allmaps(163842+1:end,s));

          % % 5 mm smoothing (sphere-based) after rank rescale
          % mapL = elec2map_apply(rankSmoothL, mapL);
          % mapR = elec2map_apply(rankSmoothR, mapR);
        case 'quantile'
          mapL = rescale_one_quantile(mapL, data.rescale_cache.L);
          mapR = rescale_one_quantile(mapR, data.rescale_cache.R);
        case 'zscore'
          mapL = rescale_one_z(mapL, data.rescale_cache.L);
          mapR = rescale_one_z(mapR, data.rescale_cache.R);
        otherwise
          error('Unknown rankrescale method: %s', cfg.rankrescale);
      end

      permMaps(:,s) = [mapL; mapR];
      if strcmp(cfg.keepsurrogates, 'yes')
        surr_subj(s,:,iter) = permMaps(:,s);
      end
    end
  else
    permMaps = allmaps; % empirical pass
  end

  % optional per-subject z-normalization before group t
  if strcmpi(cfg.normalize, 'yes')
    mu  = nanmean(permMaps, 1);
    sd  = nanstd(permMaps, 0, 1);
    permMaps = (permMaps - mu) ./ sd; % standardize
  end

  % vertex-wise one-sample t across subjects
  mu  = nanmean(permMaps, 2);
  sd  = sqrt(nanvar(permMaps, 0, 2));
  n   = double(sum(~isnan(permMaps), 2));
  tmap = mu ./ (sd ./ sqrt(n));
  tmap(isinf(tmap)) = nan;

  if strcmpi(cfg.equalize_n, 'yes')
    neff  = median(n(n>0));
    tmap  = tmap .* sqrt(neff ./ n);
  end
  tmap(mask) = nan;

  % update null or finalize observed stats
  [stat, posdistribution, negdistribution] = clusterstat(...
    tmap, sph, adj, postailcritval, negtailcritval, iter, cfg, posdistribution, negdistribution);

  if strcmp(cfg.keepsurrogates, 'yes') && iter <= cfg.numrandomization
    surr_group(:,iter) = tmap;
  end
end

% ---- Outputs
stat.coverage         = coverage;
stat.cfg              = cfg;
if strcmp(cfg.keepmaps, 'yes')
  stat.maps = allmaps;
end
if strcmp(cfg.keepsurrogates, 'yes')
  stat.surrogate_subject = surr_subj;
  stat.surrogate_group   = surr_group;
end


% =========================================================================
% Subfunctions
% =========================================================================
function [map, hemi_cache, idx, e2m_cache] = elec2map_build(side, cfg, surf_pial, surf_sph, pos_all, dv_all, sigma_mm, hemi_mask, e2m_cache)
% Per-hemisphere projector build + single apply, caching heavy aux in e2m_cache.
% side      : 'L' or 'R'
% surf_*    : hemi meshes (pial & sphere.reg)
% pos_all   : Nch x 3 electrode coords (fsavg)
% dv_all    : Nch x 1 dependent variable
% sigma_mm  : Gaussian σ in mm if cfg.kernel='gaussian'; radius(mm) if boxcar
% hemi_mask : logical Nch x 1 (channels for THIS hemi)
% e2m_cache : struct with KD-trees/aux for reuse

nVert = size(surf_pial.pos,1);
map   = nan(nVert,1);
hemi_cache = struct('P',[],'mask',[],'nVert',nVert);
idx   = int32([]);

if ~any(hemi_mask)
  return
end

% subset electrodes & DV for this hemi
posH = pos_all(hemi_mask,:);
dvH  = dv_all(hemi_mask);

% nearest pial vertices (per hemi)
idx  = int32(knnsearch(surf_pial.pos, posH, 'K',1, 'NSMethod','kdtree'));

switch lower(cfg.smooth)
  % ====================== SPHERE / EUCLIDEAN ======================
  % Fast. Builds weights using Euclidean distances in fsaverage sphere.reg
  % coordinates. Distances are in sphere units.
  case 'sphere'
    % KD-tree on sphere verts
    if side=='L'
      if isempty(e2m_cache.sphere.kdtL), e2m_cache.sphere.kdtL = createns(surf_sph.pos,'NSMethod','kdtree'); end
      kdt = e2m_cache.sphere.kdtL;
    else
      if isempty(e2m_cache.sphere.kdtR), e2m_cache.sphere.kdtR = createns(surf_sph.pos,'NSMethod','kdtree'); end
      kdt = e2m_cache.sphere.kdtR;
    end

    I=[]; J=[]; D=[];
    if strcmpi(cfg.kernel,'gaussian')
      % σ in sphere units, support ~ 3σ
      sig_unit = mm2sphere(sigma_mm);
      R        = 3*max(sig_unit, eps);  
      cand_all = rangesearch(kdt, surf_sph.pos(double(idx),:), R);
      for e=1:numel(idx)
        s    = double(idx(e));
        cand = cand_all{e}(:);  if isempty(cand), continue, end
        dvec = sqrt(sum((surf_sph.pos(cand,:) - surf_sph.pos(s,:)).^2,2));
        keep = dvec <= R; cand=cand(keep); dvec=dvec(keep);
        I=[I; cand]; J=[J; repmat(e,numel(cand),1)]; D=[D; dvec];
      end
      Wvals = exp(-(D.^2)/(2*sig_unit^2));

    elseif strcmpi(cfg.kernel,'boxcar')
      % radius in sphere units, hard cutoff
      R = max(mm2sphere(cfg.kernelwidth), eps);
      cand_all = rangesearch(kdt, surf_sph.pos(double(idx),:), R);
      for e=1:numel(idx)
        s    = double(idx(e));
        cand = cand_all{e}(:);  if isempty(cand), continue, end
        dvec = sqrt(sum((surf_sph.pos(cand,:) - surf_sph.pos(s,:)).^2,2));
        keep = dvec <= R; cand=cand(keep);
        I=[I; cand]; J=[J; repmat(e,numel(cand),1)];
      end
      Wvals = ones(size(I));
    else
      error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
    end

    % % COL
    W   = sparse(I,J,Wvals,nVert,numel(idx));
    col = full(sum(W,1));
    C   = spdiags(1./max(col,eps)', 0, size(W,2), size(W,2));
    P   = W * C; 

    % % ROW
    den = full(sum(W,2));
    m   = den>0;
    % P   = spdiags(1./max(den,eps),0,nVert,nVert) * W;

    % ====================== PIAL / GEODESIC (mm) ======================
    % Slower but more anatomically faithful. Uses geodesic distances on
    % the pial mesh (mm) via a Dijkstra traversal truncated at ~3·σ.
    % Recommended when precise folding geometry matters or kernels are
    % large. Expect ~5–10× slower projector builds on typical meshes.
  case 'pial'
    % ensure geodesic aux for this hemi is in cache
    if side=='L'
      if isempty(e2m_cache.pial.L), e2m_cache.pial.L = local_build_hemi_aux(surf_pial.pos, surf_pial.tri); end
      aux = e2m_cache.pial.L;
    else
      if isempty(e2m_cache.pial.R), e2m_cache.pial.R = local_build_hemi_aux(surf_pial.pos, surf_pial.tri); end
      aux = e2m_cache.pial.R;
    end

    I=[]; J=[]; D=[];
    if strcmpi(cfg.kernel,'gaussian')
      sig = max(sigma_mm, eps);  % σ in mm
      R   = 3*sig;
      for e=1:numel(idx)
        s = double(idx(e));
        % truncated Dijkstra to R (mm)
        d  = inf(nVert,1); d(s)=0; seen=false(nVert,1);
        heap_idx=int32(0); heap_key=0; hsz=1; heap_idx(1)=int32(s); heap_key(1)=0;
        while hsz>0
          v  = heap_idx(1); dvh = heap_key(1);
          hsz = hsz-1;
          if hsz>0
            heap_idx(1)=heap_idx(hsz+1); heap_key(1)=heap_key(hsz+1);
            p=1; while true
              l=2*p; r=l+1; m=p;
              if l<=hsz && heap_key(l)<heap_key(m), m=l; end
              if r<=hsz && heap_key(r)<heap_key(m), m=r; end
              if m==p, break; end
              [heap_idx(p),heap_idx(m)] = deal(heap_idx(m),heap_idx(p));
              [heap_key(p),heap_key(m)] = deal(heap_key(m),heap_key(p));
              p=m;
            end
          end
          if seen(v) || dvh>R, continue; end
          seen(v)=true;
          nb = aux.nbr{v}; wg = aux.wgt{v};
          for k=1:numel(nb)
            u = nb(k);
            alt = dvh + wg(k);
            if alt < d(u) && alt <= R
              d(u)=alt; hsz=hsz+1;
              if hsz>numel(heap_idx)
                newsz=max(hsz, ceil(numel(heap_idx)*1.5));
                heap_idx(end+1:newsz)=0; heap_idx=int32(heap_idx);
                heap_key(end+1:newsz)=0;
              end
              heap_idx(hsz)=u; heap_key(hsz)=alt;
              q=hsz; while q>1
                p=floor(q/2); if heap_key(p)<=heap_key(q), break; end
                [heap_idx(p),heap_idx(q)] = deal(heap_idx(q),heap_idx(p));
                [heap_key(p),heap_key(q)] = deal(heap_key(q),heap_key(p)); q=p;
              end
            end
          end
        end
        hit = find(isfinite(d) & d<=R);
        if ~isempty(hit)
          I=[I; hit]; J=[J; repmat(e,numel(hit),1)]; D=[D; d(hit)];
        end
      end
      Wvals = exp(-(D.^2)/(2*(sig^2)));

    elseif strcmpi(cfg.kernel,'boxcar')
      R = max(cfg.kernelwidth, eps);  % radius in mm; hard cutoff
      for e=1:numel(idx)
        s = double(idx(e));
        d  = inf(nVert,1); d(s)=0; seen=false(nVert,1);
        heap_idx=int32(0); heap_key=0; hsz=1; heap_idx(1)=int32(s); heap_key(1)=0;
        while hsz>0
          v  = heap_idx(1); dvh = heap_key(1);
          hsz = hsz-1;
          if hsz>0
            heap_idx(1)=heap_idx(hsz+1); heap_key(1)=heap_key(hsz+1);
            p=1; while true
              l=2*p; r=l+1; m=p;
              if l<=hsz && heap_key(l)<heap_key(m), m=l; end
              if r<=hsz && heap_key(r)<heap_key(m), m=r; end
              if m==p, break; end
              [heap_idx(p),heap_idx(m)] = deal(heap_idx(m),heap_idx(p));
              [heap_key(p),heap_key(m)] = deal(heap_key(m),heap_key(p));
              p=m;
            end
          end
          if seen(v) || dvh>R, continue; end
          seen(v)=true;
          nb = aux.nbr{v}; wg = aux.wgt{v};
          for k=1:numel(nb)
            u = nb(k); alt = dvh + wg(k);
            if alt < d(u) && alt <= R
              d(u)=alt; hsz=hsz+1;
              if hsz>numel(heap_idx)
                newsz=max(hsz, ceil(numel(heap_idx)*1.5));
                heap_idx(end+1:newsz)=0; heap_idx=int32(heap_idx);
                heap_key(end+1:newsz)=0;
              end
              heap_idx(hsz)=u; heap_key(hsz)=alt;
              q=hsz; while q>1
                p=floor(q/2); if heap_key(p)<=heap_key(q), break; end
                [heap_idx(p),heap_idx(q)]=deal(heap_idx(q),heap_idx(p));
                [heap_key(p),heap_key(q)]=deal(heap_key(q),heap_key(p)); q=p;
              end
            end
          end
        end
        hit = find(isfinite(d) & d<=R);
        if ~isempty(hit)
          I=[I; hit]; J=[J; repmat(e,numel(hit),1)];
        end
      end
      Wvals = ones(size(I));
    else
      error('cfg.kernel must be ''gaussian'' or ''boxcar''.');
    end

    W   = sparse(I,J,Wvals,nVert,numel(idx));
    den = full(sum(W,2));
    m   = den>0;
    P   = spdiags(1./max(den,eps),0,nVert,nVert) * W;

  otherwise
    error('Unknown smoothing mode "%s"', cfg.smooth);
end

% apply once, return cache for reuse
tmp = P * dvH; map(m) = tmp(m);
hemi_cache.P    = P;
hemi_cache.mask = m;

% -------------------------------------------------------------------------
function map = elec2map_apply(hemi_cache, dv_hemi)
% Reuse projector/mask for permutations.
map = nan(hemi_cache.nVert,1);
if ~isempty(hemi_cache.P)
  tmp = hemi_cache.P * dv_hemi;
  if ~isempty(hemi_cache.mask)
    map(hemi_cache.mask) = tmp(hemi_cache.mask);
  else
    map = tmp;
  end
end

% -------------------------------------------------------------------------
function aux = local_build_hemi_aux(V, F)
nV = size(V,1);
E  = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1]); F(:,[2 1]); F(:,[3 2]); F(:,[1 3])];
dV = V(E(:,1),:) - V(E(:,2),:);
W  = sqrt(sum(dV.^2,2));
A  = sparse(E(:,1),E(:,2),W,nV,nV); A = min(A, A.');
[i,j,w] = find(A);
nbr = accumarray(i, j, [nV 1], @(x){int32(x)});
wgt = accumarray(i, w, [nV 1], @(x){x});
aux = struct('nbr',{nbr}, 'wgt',{wgt});

% -------------------------------------------------------------------------
function su = mm2sphere(mm)
% Approximate conversion: 1 mm (pial) ≈ 0.5 units (sphere.reg)
su = 0.5 * mm;

% -------------------------------------------------------------------------
function [Q, lam] = graph_eigenbasis(pos, sigma_graph)
% (a) Build Gaussian-weighted graph (no self-loops), normalized Laplacian
% (b) Return eigenvectors (columns) ordered by eigenvalue + eigenvalues
%
% sigma_graph (mm) controls the width of the electrode-level Gaussian
% weighting on the subject’s *electrode graph* (used to form a normalized
% Laplacian and its eigenbasis for spectral surrogates). We first convert
% this mm scale to sphere units per hemi (same local scale factor as
% above), then build W_ij = exp(-d_ij^2 / (2*σ_graph^2)) among electrodes.
% Larger sigma_graph → more low-frequency modes (smoother surrogates).
n = size(pos,1);
if n==0
  Q = zeros(0,0); lam = []; return
elseif n==1
  Q = 1; lam = 0; return
end
r = sqrt(sum(pos.^2, 2));
U = pos ./ r;  % unit vectors
C = U * U.';                 
C = max(-1, min(1, C));
theta = acos(C);             % radians
D = 100 * theta;
W  = exp(-(D.^2) ./ (2*(sigma_graph^2)));
W(1:n+1:end) = 0; % prevent self-loops
d  = sum(W,2); d(d==0) = 1;
Dm = spdiags(1./sqrt(d), 0, n, n);
L  = speye(n) - (Dm * W * Dm);
L  = (L + L')/2;
[Qtmp, Lam] = eig(L);
lam = diag(Lam);
[lam, ord] = sort(lam, 'ascend'); % smooth → rough
Q = Qtmp(:,ord);

% -------------------------------------------------------------------------
function a2 = spectral_bandrotate(a, lam)
% Random orthogonal rotation within eigenvalue bands (preserves band power)
% Inputs:
%   a   : modal coefficients (length m)
%   lam : eigenvalues aligned to a (length m)
% Output:
%   a2  : rotated coefficients (same length)
if isempty(a) || isempty(lam), a2=a; return, end
B=20; edges=quantile(lam, linspace(0,1,B+1)); a2=a;
for b=1:B
  if b<B, idx = find(lam>=edges(b) & lam<edges(b+1));
  else,    idx = find(lam>=edges(b) & lam<=edges(b+1));
  end
  k=numel(idx);
  if k>1
    G=randn(k,k); [Q,~]=qr(G,0); if det(Q)<0, Q(:,1)=-Q(:,1); end
    a2(idx)=Q*a(idx);
  end
end

% -------------------------------------------------------------------------
function out = rescale_one_quantile(in, c)
out=in; m=c.mask;
if any(m)
  vals=in(m);
  [~,ord]=sort(vals,'ascend');
  tmp=vals;
  tmp(ord)=c.tgt_sorted;
  out(m)=tmp;
end

% -------------------------------------------------------------------------
function out = rescale_one_z(in, c)
out=in; m=c.mask;
if any(m)
  x=in(m);
  x=(x-mean(x))./max(std(x),eps);
  out(m)=c.mu + c.sd * x;
end

% -------------------------------------------------------------------------
function [stat, posdistribution, negdistribution] = clusterstat(tmap, sph, adj, postailcritval, negtailcritval, iter, cfg, posdistribution, negdistribution)
if iter <= cfg.numrandomization
  % accumulate null distributions
  posclusrnd = clusterfunc(tmap >= postailcritval, sph, adj);
  negclusrnd = clusterfunc(tmap <= negtailcritval, sph, adj);
  if any(posclusrnd(:))
    poscstats = accumarray(posclusrnd(posclusrnd>0), tmap(posclusrnd>0));
    posdistribution(iter) = max(poscstats);
  end
  if any(negclusrnd(:))
    negcstats = accumarray(negclusrnd(negclusrnd>0), tmap(negclusrnd>0));
    negdistribution(iter) = min(negcstats);
  end
  stat = [];
else
  % finalize observed
  prb_pos = ones(size(tmap)); prb_neg = ones(size(tmap));

  % positive clusters
  posclusobs = clusterfunc(tmap >= postailcritval, sph, adj);
  if any(posclusobs(:))
    poscstats = accumarray(posclusobs(posclusobs>0), tmap(posclusobs>0));
    [sorted_cstats, sort_idx] = sort(poscstats, 'descend');
    map = zeros(max(posclusobs(:)),1); map(sort_idx) = 1:numel(sorted_cstats);
    posremap = nan(size(tmap)); valid = posclusobs > 0;
    posremap(valid) = map(posclusobs(valid)); posclusobs = posremap;
    for j = 1:numel(sorted_cstats)
      posclusters(j).clusterstat = sorted_cstats(j);
      posclusters(j).prob = (sum(posdistribution >= sorted_cstats(j)) + 1) / (numel(posdistribution) + 1);
    end
    prb_pos(valid) = [posclusters(posclusobs(valid)).prob];
  else
    posclusters = [];
  end

  % negative clusters
  negclusobs = clusterfunc(tmap <= negtailcritval, sph, adj);
  if any(negclusobs(:))
    negcstats = accumarray(negclusobs(negclusobs>0), tmap(negclusobs>0));
    [sorted_cstats, sort_idx] = sort(negcstats, 'ascend');
    map = zeros(max(negclusobs(:)),1); map(sort_idx) = 1:numel(sorted_cstats);
    negremap = nan(size(tmap)); valid = negclusobs > 0;
    negremap(valid) = map(negclusobs(valid)); negclusobs = negremap;
    for j = 1:numel(sorted_cstats)
      negclusters(j).clusterstat = sorted_cstats(j);
      negclusters(j).prob = (sum(negdistribution <= sorted_cstats(j)) + 1) / (numel(negdistribution) + 1);
    end
    prb_neg(valid) = [negclusters(negclusobs(valid)).prob];
  else
    negclusters = [];
  end

  % two-tailed (separate nulls)
  stat = struct();
  stat.stat                = tmap;
  stat.prob                = min(prb_pos, prb_neg); % x2 for Bonferroni correction
  stat.prob(stat.prob > 1) = 1;
  stat.mask                = stat.prob <= 0.025;  % assuming 2-tailed w/ alpha = 0.05

  stat.posclusters         = posclusters;
  stat.posclusterslabelmat = posclusobs;
  stat.posdistribution     = posdistribution;
  stat.negclusters         = negclusters;
  stat.negclusterslabelmat = negclusobs;
  stat.negdistribution     = negdistribution;
end

% -------------------------------------------------------------------------
function clustermap = clusterfunc(mask, ctx, adjacency)
sig_idx = find(mask);
if isempty(sig_idx)
  clustermap = nan(size(ctx.pos,1),1);
  return;
end
subgraph = adjacency(sig_idx, sig_idx);
G = graph(subgraph);
bins = conncomp(G);
clustermap = nan(size(ctx.pos,1),1);
clustermap(sig_idx) = bins;
