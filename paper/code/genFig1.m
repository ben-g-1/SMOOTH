function genFig1
% ------------------------------------------------------------
% Figure 1. Intracranial EEG Task and Representative Findings
%   1) Electrode-level statistical maps
%   2) ROI-level significance maps
%   3) Group-level coverage map
%   4) Exemplar channel time series
% ------------------------------------------------------------

%% config

clear; clc; close all;

% set paths
paths.smooth      =  'C:\Users\bgrau\GitHub\SMOOTH';        % UPDATE
paths.fieldtrip   = 'C:\Matlab\toolboxes\fieldtrip';     % UPDATE
paths.freesurfer  = 'C:\Matlab\toolboxes\fieldtrip\external\freesurfer';    % UPDATE
paths.demoSource  = fullfile(paths.smooth, 'demo', 'source.mat');
paths.figOut      = fullfile(paths.smooth, 'paper', 'figs', '1');
paths.tsFile      = fullfile(paths.smooth, 'paper', 'ts.mat');

% add paths
addpath(genpath(paths.smooth));
addpath(paths.fieldtrip);
ft_defaults;

% cortical mesh (only for visualization)
cortex = load_fsaverage_mesh(paths.freesurfer);

% read in data
load(paths.demoSource)  % path to example dataset
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% organize channel-level data

CLstats = struct('sub',[],'stat',[],'fsxyz',[],'fssurfxyz',[],'dk',{{}});

for i = 1:numel(source)
    CLstats.sub       = vertcat(CLstats.sub,i.*ones(length(source{i}.label),1));
    CLstats.stat      = vertcat(CLstats.stat,source{i}.stat);
    CLstats.fsxyz     = vertcat(CLstats.fsxyz,source{i}.elec.chanpos);
    
    % correct jitter
    D = pdist2(cortex.pos,source{i}.elec.chanpos);
    [~, idx] = min(D,[],1);
    fssurfxyz = cortex.pos(idx,:);
    CLstats.fssurfxyz = vertcat(CLstats.fssurfxyz,fssurfxyz);

    % ID channels
    [roi_name, ~, ~, ~] = fsavg_coord_to_dk(fssurfxyz, paths.freesurfer);
    CLstats.dk        = vertcat(CLstats.dk,roi_name);
end


%% visualize channel-level stats

% right-tail significance
ppos = 1 - normcdf(CLstats.stat);
[sig, ~, ~, ~]=fdr_bh(ppos);  % FDR correct

% colormap (orange/neutral/blue)
blue     = slanCM('Blues');
orange   = slanCM('Oranges');
berrymap = [orange(180,:); 0 0 0; blue(180,:)];

% plot
fig = figure('units','normalized','outerposition',[0 0 1 1],'visible','on','color',[1 1 1]); hold on; 
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
ft_plot_cloud(CLstats.fssurfxyz((sig & CLstats.stat>0),:), 1.*ones(sum(sig & CLstats.stat>0),1), 'cloudtype', 'surf', ...
  'radius', 2.5, 'scalerad', 'no', 'colormap', berrymap, 'clim', [-1 1]);   % sig
ft_plot_cloud(CLstats.fssurfxyz((~sig | CLstats.stat<0),:), zeros(sum(~sig | CLstats.stat<0),1), 'cloudtype', 'surf', ...
  'radius', 1, 'scalerad', 'no', 'colormap', berrymap, 'clim', [-1 1]);   % null

% export
% views    = [90 0; -90 0; 0 90; 180 -90; 180 -15; 0 0];
% angle    = {'right', 'left', 'dorsal', 'ventral', 'rostral', 'caudal'};
% for v = 1:size(views,1)
%     view(views(v,:))
%     exportgraphics(fig, [paths.figOut filesep 'berries_' angle{v} '.pdf'], 'contenttype', 'image', 'resolution', 450)
% end
% close(fig)

%% run ROI-based stats

rois  = unique(CLstats.dk);
tvals = []; 
pvals = []; 

for r = 1:length(rois)
    roi    = rois{r};
    in_roi = strcmp(roi,CLstats.dk);
    zchans = CLstats.stat(in_roi);

    sidx   = CLstats.sub(in_roi);
    subs   = unique(sidx);
    nsub   = length(subs);

    if nsub > 3
        % average w/in subs
        zsubs  = zeros(nsub,1);
        for s = 1:nsub
            zsubs(s) = mean(zchans(sidx==subs(s)));
        end
        [~,p,~,stat] = ttest(zsubs);
        tvals(end+1) = stat.tstat;
        pvals(end+1) = p;
    else
        tvals(end+1) = nan;
        pvals(end+1) = nan;
    end
end
[~,~,~,pvals] = fdr_bh(pvals);

% visualize
plot_roi_dk('rh', rois, tvals, pvals, paths.freesurfer, blue);
plot_roi_dk('lh', rois, tvals, pvals, paths.freesurfer, blue);

%% coverage map

% smooth
cfg                  = [];
cfg.fshome           = '/Applications/freesurfer/7.3.2';
cfg.numrandomization = 0;
cfg.kernelwidth      = 30;
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
stat = SMOOTHstat(cfg, source{:});  
y = stat.coverage;
y(y<3) = nan;  % coverage threshold

% visualize
f = figure('color',[1 1 1]);
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', y);
clim([0 20]); colormap(slanCM('Greens')); view([90 0])
% exportgraphics(f, [paths.figOut filesep 'cvgRH.pdf'], 'contenttype', 'image', 'resolution', 450)
% view([-90 0])
% exportgraphics(f, [paths.figOut filesep 'cvgLH.pdf'], 'contenttype', 'image', 'resolution', 450)
% close(f)


%% plot colorbars
f = figure;
axis off
hcb = colorbar;
clim([0 20])
hcb.Ticks = 0:20;
colormap(slanCM('Greens'))
exportgraphics(f, [paths.figOut filesep 'greenCB.pdf'], 'contenttype', 'vector')
colormap(cmap(1:end-1,:)./255)
hcb.Ticks = [];
% exportgraphics(f, [paths.figOut filesep 'blorangeCB.pdf'], 'contenttype', 'vector')
% close(f)

%% plot colorbars
f = figure;
axis off
hcb = colorbar;
hcb.Ticks = 0:5;
clim([0 5])
colormap(slanCM('Blues'))
% exportgraphics(f, [paths.figOut filesep 'blueCB.pdf'], 'contenttype', 'vector')

%% exemplar channel timeseries
load([paths.smooth filesep 'paper/ts.mat'])
f = figure('color',[1 1 1],'Position',[100 300 1500 415]);
plot(ts.t,ts.y);
hold on; plot(ts.t,ts.r);
xlim([0 300])
xticks(0:150:300)
ylim([-1 2])
yticks(-1:2)
xlabel('time (s)')
ylabel('zHGB')
legend({'zHGB','mentalizing'})
% exportgraphics(f, [paths.figOut filesep 'timeseries.pdf'], 'contenttype', 'vector')
% close(f)

end

%% subfunctions
% -------------------------------------------------------------------------
function cortex = load_fsaverage_mesh(fshome)
% load cortex mesh given fs path

lh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','lh.pial'));
rh = ft_read_headshape(fullfile(fshome,'subjects','fsaverage','surf','rh.pial'));

cortex      = [];
cortex.pos  = [lh.pos; rh.pos];
cortex.tri  = [lh.tri; rh.tri + size(lh.pos,1)];
cortex.unit = lh.unit;

if isfield(lh,'curv') && isfield(rh,'curv')
    cortex.curv = [lh.curv; rh.curv];
else
    cortex.curv = zeros(size(cortex.pos,1),1);
end
end
% -------------------------------------------------------------------------
function idx = norm_to_cmap_idx(x, xmin, xmax, nColors)
%NORM_TO_CMAP_IDX  Map values to 1..nColors; NaNs -> nColors+1.
idx = floor((x - xmin) ./ (xmax - xmin) * nColors) + 1;
idx(idx < 1) = 1;
idx(idx > nColors) = nColors;
idx(isnan(idx)) = nColors + 1;
end
% -------------------------------------------------------------------------
function plot_roi_dk(hemi, rois, tvals, pvals, fshome, blue)
% render ROI-level -log(p) on DK parcellation (one hemisphere).
assert(ismember(hemi, {'lh','rh'}), 'hemi must be ''lh'' or ''rh''.');

annotFile = fullfile(fshome,'subjects','fsaverage','label',sprintf('%s.aparc.annot', hemi));
[~,~,actbl] = read_annotation(annotFile);
dk_order = actbl.struct_names;

isHemi = cellfun(@(x) contains(x,hemi), rois);
hemi_rois = rois(isHemi);
hemi_t    = tvals(isHemi);
hemi_p    = pvals(isHemi);

% strip 'lh_'/'rh_' prefix to match annotation struct_names
hemi_rois = cellfun(@(x) x(4:end), hemi_rois, 'UniformOutput', false);

% map into DK order (prepend NaN for 'unknown')
hemi_p = [NaN; hemi_p(:)];
hemi_t = [NaN; hemi_t(:)];

[~, order] = ismember(dk_order, hemi_rois);
plotp = -log(hemi_p(order+1));
plott = hemi_t(order+1);

% normalize -log(p) into 1..256, NaNs -> 257 (background)
idx = norm_to_cmap_idx(plotp, 0, 5, 256);
idx(plott < 0) = 257;

bg   = [227/255 222/255 193/255];
cmap = 255 .* [blue; bg];  % 256 blues + background
cfg = [];
cfg.elecCoord = 'n';
cfg.surfType  = 'pial';
cfg.overlayParcellation = 'DK';
cfg.fsurfSubDir = fullfile(fshome,'subjects');
cfg.view = hemi(1); % 'l' or 'r'
cfg.parcellationColors = cmap(idx, :);

plotPialSurf('fsaverage', cfg);
end
% -------------------------------------------------------------------------
