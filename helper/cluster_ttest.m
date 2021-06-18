function [stat, ftdata] = cluster_ttest(data, srate, t_range, neighbours, ...
    labchans, alpha, t_calc, eeg_chan, num_tails)
%
% wrapper function for doing cluster-based erp analysis using field trip
% the code is adapted from fieldtriptoolbox.org/tutorial/cluster_permutation_timelock
%
% data: two-cell containing arrays for doing ttest where first dimension is subject
% followed by channel x time-frequency
% srate: sampling rate
% t_range: specifies the time values for the erp
% neighbiours: neighbour structure
% labchans: channel labels to analyse
% t_calc: time range over which to calculate the clustering

if ~exist('t_calc', 'var'), t_calc = [0 t_range(end)]; end
if ~exist('eeg_chan', 'var'), eeg_chan = true; end
if ~exist('alpha', 'var')||isempty(alpha), alpha = .025; end
if ~exist('num_tails', 'var')||isempty(num_tails), num_tails = 2; end

nsubj = size(data{1},1);
allsubjyes = cell(1, nsubj);
allsubjno = allsubjyes;

for ns = 1:nsubj
    allsubjyes{ns}.label = labchans';
    allsubjyes{ns}.fsample = srate;
    allsubjyes{ns}.avg = permute(data{1}(ns, :, :, :, :, :), [2:6 1]);
    
    allsubjyes{ns}.time = t_range;
    if eeg_chan
        allsubjyes{ns}.dimord = 'chan_time';
    else
        allsubjyes{ns}.dimord = 'time';
    end
    
    allsubjno{ns}.label = labchans';
    allsubjno{ns}.fsample = srate;
    allsubjno{ns}.avg = permute(data{2}(ns, :, :, :, :, :), [2:6 1]);
    
    allsubjno{ns}.time = t_range;
    if eeg_chan
        allsubjno{ns}.dimord = 'chan_time';
    else
        allsubjno{ns}.dimord = 'time';
    end
end

%

cfg         = [];

if eeg_chan
    cfg.channel = {'EEG'};
end
cfg.latency = t_calc;

cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 1;
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = alpha;
cfg.numrandomization = 'all';

design = zeros(2, nsubj*2);
design(1,:) = [1:nsubj 1:nsubj];
design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

[stat] = ft_timelockstatistics(cfg, allsubjyes{:}, allsubjno{:});

if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    pos_cluster_pvals = [stat.posclusters(:).prob];
    pos_clust = find(pos_cluster_pvals < .05/num_tails);
    pos       = ismember(stat.posclusterslabelmat, pos_clust);
    figure, subplot(211), imagesc(pos), title('positive clusters')
end

if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_clust = find(neg_cluster_pvals < .05/num_tails);
    neg       = ismember(stat.negclusterslabelmat, neg_clust);
    subplot(212), imagesc(neg), title('negative clusters')
end
ftdata{1} = allsubjyes;
ftdata{2} = allsubjno;

end