% -------------------------------------------------------------------------
% demo script for SMOOTH group-level stats
%
% demonstrates a minimal end-to-end workflow:
%   1) set up paths (toolbox + FieldTrip)
%   2) load an example FieldTrip "source" dataset
%   3) run SMOOTHstat (group-level inference)
%   4) visualize the resulting t-map and significance mask
%   5) evaluate surrogate data
%
% requirements:
%   - FieldTrip
%   - FreeSurfer installation (for fsaverage surfaces)
%   - example dataset: demo/source.mat (distributed with the toolbox)
%
% Authors: Kamren Khan & Arjen Stolk, v0.1
% -------------------------------------------------------------------------

clear; clc; close all;

% set paths
paths.smooth      =  'C:\Users\bgrau\GitHub\SMOOTH';        % UPDATE
paths.fieldtrip   = 'C:\Matlab\toolboxes\fieldtrip';     % UPDATE
paths.freesurfer  = 'C:\Matlab\toolboxes\fieldtrip\external\freesurfer';    % UPDATE

% add paths
addpath(genpath(paths.smooth))
addpath(paths.fieldtrip)  
ft_defaults

% cortical mesh (only for visualization)
cortex      = load_fsaverage_mesh(paths.freesurfer);

disp(paths.smooth)
% read in data
load([paths.smooth filesep 'demo' filesep 'source.mat'])  % path to example dataset  % making both filesep helps for Windows. BG
nsubs = numel(source);
nchans = sum(cellfun(@(s) numel(s.label), source));
fprintf('\nLoaded demo dataset: %d subjects, %d total electrodes.\n', nsubs, nchans);

%% run SMOOTH
tic
cfg                  = [];
cfg.keepmaps         = 'yes';
cfg.keepsurrogates   = 'yes';
cfg.normalize        = 'no';
cfg.numrandomization = 500;  % created array with size above 32 GB at 500 randomizations
cfg.kernelwidth      = 35;
cfg.rankrescale      = 'quantile';
cfg.smooth           = 'sphere';
cfg.kernel           = 'gaussian';
cfg.fshome           = paths.freesurfer;
stat = SMOOTHstat(cfg, source{:}); %threw error: 
                                        % "Unable to perform assignment because the left and right sides have a different number of elements.

                                        % Error in SMOOTHstat (line 210)
                                        %           [~, idxL] = sort(mapL); mapL(idxL) = sort(allmaps(1:163842,s));

                                        % Error in SMOOTHdemo (line 53)
                                        % stat = SMOOTHstat(cfg, source{:});"
t = toc;
fprintf('SMOOTHstat completed in %.2f s.\n', t);

%% visualize output

% (1) group-level t-map
figure
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', stat.stat);
clim([-3 3]); colormap('parula')
title('group-level t-map'); view([90 0])

% (2) significant effects
figure
y = stat.stat;  % pos effects
y(~stat.mask | stat.stat<0) = nan;
ft_plot_mesh(cortex, 'vertexcolor', 'curv');
hold on; ft_plot_mesh(cortex, 'vertexcolor', y);
clim([-3 3]); colormap('parula')
title('SMOOTH significant effects'); view([90 0])

%% surrogate robustness check
% verify that surrogate spatial autocorr is equivalent to that of empirical
% data

% surrogate diagnostics
d    = SMOOTHdiag(stat); % Lack of verbosity made me assume this line hung. Why are all the fprintf lines in the function commented out? BG
emp  = d.subject.spatial.BOTH.moranI.per_subject.emp(:);
null = d.subject.spatial.BOTH.moranI.per_subject.null_mu(:);
[~,p,~,tstat] = ttest(emp, null);  % emp vs null smoothness

% vis
figure; hold on;
scatter(ones(size(emp)), emp, 25, 'filled');
scatter(2*ones(size(null)), null, 25, 'filled');
plot([ones(numel(emp),1) 2*ones(numel(emp),1)]', [emp null]', 'k-');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Empirical','Null'});
title("spatial autocorrelation robustness check");
subtitle(sprintf('paired t-test: t = %.2f, p = %.3g', tstat.tstat, p));
box off;
