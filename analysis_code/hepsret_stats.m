function hepsret_stats_hep(exp_hepsret)

% load 32 channel locations
cd(fullfile(exp_hepsret.session_dir))

load('chanlocs32.mat')
load('hrep.mat')

% compute neighbour structure
neighbourdist = 50;
neighbours = neighbours_struct(chanlocs, neighbourdist);

%%% analyse hep 
% calculate cluster interaction effect
data{1} = squeeze(diff(hrep(:, 1, :, :, :), [], 3));
data{2} = squeeze(diff(hrep(:, 2, :, :, :), [], 3));
nht = size(hrep,5);
tt = linspace(trange_hep(1), trange_hep(2), nht);
[stat, ftdata] = cluster_ttest(data, srate, tt, neighbours, labchans);

% plot interaction effect
timestep      = 0.02; %(in seconds)
cluster_tplot(ftdata, stat, srate, tt, timestep)

% for the interaction clusters do t-test separately for positive and
% negative channels
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);
figure, imagesc(pos), title('positive clusters')

h2 = hrep(:,:,:,:,tt>=0);
tpos = logical(sum(pos,1));
cpos = logical(sum(pos,2));
h2 = mean(mean(h2(:,:,:,cpos,tpos),4),5);

[h,p] = ttest(squeeze(diff(h2(:,1,:),[],3))); % pos words 
[~,p,~,stat] = ttest(squeeze(diff(h2(:,2,[2 1]),[],3))); % neg words
fprintf('neg words yes > no, t(%g)=%g, p=%g\n', stat.df, stat.tstat, p)


%%%
%%% calculate cluster main effect
data{1} = squeeze(sum(hrep(:, :, 1, :, :), 2));
data{2} = squeeze(sum(hrep(:, :, 2, :, :), 2));
[stat, ftdata] = cluster_ttest(data, srate, tt, neighbours, labchans);


end