%% Data analysis script to investigate
%% whether prestimulus pulse-evoked potentials (PEPs) influence
%% self-association with affective (pleasant and unpleasant) adjectives
%%

%% The dataset has been previously published with regards to a different
%% question: Katyal et al (2020), CABN
%%


%% add necessary paths

load_hepsret
add_fieldtrip_path
rmpath ~/OneDrive/Projects/NeuroMatlabToolboxes/fieldtrip-20210310/external/signal/


%% Initialize experiment parameters

exp_hepsret = hepsret_experiment;
filename = 'csd_icacomprem_sret4_medriv_physio';

%%
%%
%%
%% Step 1.1 - Get the average reaction times from the SRET experiment
%%% based on this decide the different values of the intervals over which
%%% to calculate the PEPs
hepsret_rt(exp_hepsret)

%%
%% Step 1.2 - Perform simulations 
%%% to determine if the dataset has enough power to test the hypotheses.
%%% This is needed because there was a self-positivity bias in responses,
%%% leaving us with far fewer trials for the pleasant-nonself vs.
%%% pleasant-self and unpleasant-self vs. unpleasant-nonself conditions

% iterate over the different values of the interval over which to determine
% pulse peaks for calculating PEPs
bef_hep_iter = -3:-1:-6;
nbef = numel(bef_hep_iter);


% iterate over different values of the positivity-bias criterion. a minimum
% of 2,3,4 and 5 trials for pleasant-nonself and unpleasant-self conditions 
endorscutoff_range = [.84 .87 .90 .94 .97];
necr = numel(endorscutoff_range);

 % duration of the PEP, baseline is subtracted between first index and time
 % zero
trange_hep = [-.1 .5];

% estimate the effect size based on the difference between high- and
% low-arousal conditions in Luft and Bhattacharya (2016) Sci Rep
% Figure 3 (in µV/m3) as well as Marshall et al (2018) SCAN
effect_size = 3.75;

% number of permutations
nperms = 10000; 
alpha = 0.05;

% initialse simulated power
probdetect = NaN(necr, nbef, 2);
% initialise number of included subjects for each bias criterion
nsubj_valid = NaN(necr, nbef);


% loop over different interval durations
for nb = 1:nbef
    bef_aft_hep = [bef_hep_iter(nb) 0];
    
    
    % loop over the different bias criteria
    for ecr = 1:necr
        endorscutoff = endorscutoff_range(ecr);
        
        % first determine the mean number of PEP trials per subject
        % and per Valence (pleasant, unpleasant) x Endorse (Yes, No)
        % conditions fromt the dataset itself. also determine the average
        % noise in EEG 
        % channels (defined as stdev of frequencies > 7 Hz) across
        % subjects, and the number of subjects for each threshold
        [hbtrials, noise_subj, val_subj] = hepsret_optim_inclus...
            (exp_hepsret, filename, ...
            bef_aft_hep, trange_hep, endorscutoff);
        nsubj_valid(ecr,nb) = sum(val_subj);
        
        % loop over pleasant and unpleasant words
        for pn = 1:2 
            % simulate the experiment with expected effect size and
            % average noise for each threshold
            probdetect(ecr,nb,pn) = hepsret_simulate_detect_prob...
                (effect_size, noise_subj, squeeze(round(hbtrials(:, pn, :))),...
                nsubj_valid(ecr,nb), nperms, alpha);
        end
    end
end

cd(exp_hepsret.session_dir)
save('hepsret_datafiles/simulate_power', 'nbef', 'probdetect',...
    'nsubj_valid', 'bef_hep_iter', 'nperms', 'effect_size')

%% Step 1.3 - plot simulations

load('hepsret_datafiles/simulate_power')
figure(1)
for nb = 1:nbef
    subplot(2,2,nb)
    xx = 100-(5:-1:1)/30*100;
    plot(xx, squeeze(probdetect(:,nb,:)), '-.', 'LineWidth', 1.5)
    if nb==1
        text(xx+[.9 1 -2.1 .9 -2.6], ...
            probdetect(:,nb,1)'+[.002 0 .005 -.02 0],...
            {'n = 6', 'n = 9', 'n = 12', 'n = 13', 'n = 19'})
    end
    hold on,plot(xx,mean(probdetect(:,nb,:),3), '-ok', 'LineWidth', 4)
    title(sprintf('interval = %g sec', bef_hep_iter(nb)))
    ax = [82 98 0.4 .95]; axis(ax)
    legend({'Pleasant', 'Unpleasant', 'Average'})
end

%% Step 1.4 - Selected parameters: interval = 5 sec, bias criterion = 90%

%%% recalculate the valid participants for the selected values
endorscutoff = 0.9;
bef_aft_hep = [-5 0];

[hbtrials, ~, val_subj] = hepsret_optim_inclus...
    (exp_hepsret, filename, ...
    bef_aft_hep, trange_hep, endorscutoff);
selsubj = exp_hepsret.data(val_subj);
subjage = [selsubj.subj_age];
fprintf('age of selected participants: mean=%g, stdev=%g\n', ...
    mean(subjage), std(subjage))
fprintf('number of females: %g\n', sum([selsubj.subj_sex]))

hbt2 = reshape(hbtrials, [12 4]);
mean(hbt2,1)
std(hbt2,1)

%%
%%
%%
%% Step 2 - Main analysis step 

%%% Parse through the dataset and get the PEPs, ERPs, and HR-ERPs and save
%%% them in a separate file

selchans = {};
 
bef_aft_hep = [-5 0]; % interval over which to evaluate the HEP
trange_hep = [-.1 .5]; % HEP window
bef_aft_erp = [-.2 .6]; % ERP window
bef_aft_hr = [-5 0]; % window for calculating HR

smooth_erp = 0; % in ms 
blsub_flag = true; % subtract baseline for ERP and HEP
avg_chans = false; % don't average channels, return them separately

endorscutoff = .9; %%% winning bias criterion based on simluations

shuffle_peaks = 0;

%%% to plot the detected peaks for each trial (green and red circles for
%%% included and excluded peaks respectively set do_plots = 1
do_plots = 0; 

hepsret_hep(exp_hepsret, selchans, filename, ...
    smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr,...
    blsub_flag, endorscutoff, shuffle_peaks, do_plots);


%%
%%
%%
%% Step 3 - Stats

%% 3.0 - Compute the neighbours structure for performing cluster-based analysi

cd(exp_hepsret.session_dir)
load('chanlocs.mat')
neighbourdist = 50;
neighbours = neighbours_struct(chanlocs, neighbourdist);

%% 3.1.1 Cluster-based stats - Pleasant word PEP

stat_type = 2; % for plotting pleasant HEP
hepsret_stats(exp_hepsret, stat_type, neighbours);

%%% no significant clusters for pleasant word heps for yes > no

%% 3.1.2 Cluster-based stats - Unpleasant word PEP

stat_type = 3; % for plotting unpleasant HEP
[stats_hep, ftdata] = hepsret_stats(exp_hepsret, stat_type, neighbours);
fprintf('p value of unpleasant word effect: %g\n', stats.posclusters(1).prob)

%% 3.1.3 Plot topography of the effect in 3.1.2

timestep = 0.024; % in seconds
srate = 250;
plot_starstop = [12 13]; % plot topography of only a subset of time points
cluster_tplot(ftdata, stats_hep, srate, timestep, [])

%% 3.1.4 Plot HEPs across all channels

plot_type = 1; % 1-hep, 2-erp, 3-hr
hepsret_ploterps(exp_hepsret, plot_type)

%% 3.1.5 Plot Figure 3 with significant channels HEPs and topography

fig_num = 1; % 
hepsret_figures(exp_hepsret, stats_hep, fig_num)

%%
%% Control analyses
%%
%% 3.2.0 Cluster-based stats - Unpleasant word Heartrate (control analysis)

stat_type = 7; % cluster stats on heart-rate
hepsret_stats(exp_hepsret, stat_type, neighbours);

% found no clusters for difference in prestimulus heart rate for unpleasant
% words 

%% 3.2.1 - Shuffling

shuffle = 2000; % number of shuffles
hepsret_hep(exp_hepsret, selchans, filename, ...
    smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr,...
    blsub_flag, endorscutoff, shuffle);

%% 3.2.2 Plot shuffled data at significant channels from main analysis of Unpleasant words

fig_num = 2; % 
hepsret_figures(exp_hepsret, stats_hep, fig_num)

%% 3.2.3 Do a t-test of shuffled at the significant channel-by-time cluster from main analysis

stat_type = 2; % HEPs for unpleasant words
reg_or_shuf = 2; % 1-regular unshuffled data, 2-shuffled data
% this function takes in a stat
hepsret_t_hepdur(exp_hepsret, stats_hep, stat_type, reg_or_shuf);

%% 3.2.4 Cluster-based stats - Unpleasant word HR - Shuffled

stat_type = 12; % for stats on shuffled negative HEP
[stats_shuf, ftdata] = hepsret_stats(exp_hepsret, stat_type, neighbours);

%% try to do logistic regression on individual trials within subjects

% hepsret_logistic_regression(exp_hepsret, stats, filename,...
%     bef_aft_hep, trange_hep, bef_aft_hr,...
%     blsub_flag, endorscutoff)
% 
% %%% doesn't work because of the extremity in bias

%     %% do t-ttest of the cluster from stats to the currently analysed hrep.mat file
%     stat_type = 2; % 1-pos HEP, 2-neg HEP, 3-pos HR, 4-neg HR
%     hepsret_t_hepdur(exp_hepsret, stats, stat_type)
%
%     %% do t-ttest of the cluster from stats to the currently analysed hrep.mat file
%     stat_type = 4; % 1-pos HEP, 2-neg HEP, 3-pos HR, 4-neg HR
%     hepsret_t_hepdur(exp_hepsret, stats, stat_type)

%%
%%
%% 3.3.1 How does the significant effect change with prestimulus interval duration
%%%  to a t-test for negative words HEP with different values of the HEP
%%%  evaluation interval

interv_start_dur = -10:.25:-1.5;
nintervs = numel(interv_start_dur);
stat_type = 2; % 1-pos HEP, 2-neg HEP, 3-pos HR, 4-neg HR
tvals = NaN(1, nintervs); pvals = tvals;

for nin = 1:nintervs
    bef_aft_hep = [interv_start_dur(nin) 0];
    hepsret_hep(exp_hepsret, selchans, filename, ...
        smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr,...
        blsub_flag, endorscutoff,...
        do_plots);
    [tvals(nin), pvals(nin)] = hepsret_t_hepdur(exp_hepsret, stats, stat_type);
end
save('hepsret_datafiles/t_diff_intervals', 'interv_start_dur', 'tvals', ...
    'pvals')

%% 3.3.2 - Plot Figure 4a - t vs. interval duration

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/t_diff_intervals')
figure,plot(-interv_start_dur,tvals, 'LineWidth', 3)
hold on, plot([-interv_start_dur(1) 1], [2.2 2.2], 'k-.', ...
    'LineWidth', 2)

%% 3.3.2 How does the significant effect change with prestimulus interval centre
%%% t-test for negative words HEP with different centres of the HEP
%%% evaluation interval 
    
interv_mid_dur = -4.5:.5:-1.5;
nintervs = numel(interv_mid_dur);
interv_dur = 3;
stat_type = 2; % 1-pos HEP, 2-neg HEP, 3-pos HR, 4-neg HR
tvals_mid = NaN(1, nintervs); pvals_mid = tvals_mid; dvals_mid = NaN(12, nintervs);

for nin = 1:nintervs
    bef_aft_hep = [interv_mid_dur(nin)-interv_dur/2 ...
        interv_mid_dur(nin)+interv_dur/2];
    hepsret_hep(exp_hepsret, selchans, filename, ...
        smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr,...
        blsub_flag, endorscutoff,...
        do_plots);
    [tvals_mid(nin), pvals_mid(nin), dvals_mid(:,nin)] = ...
        hepsret_t_hepdur(exp_hepsret, stats, stat_type);
end

save('hepsret_datafiles/t_diff_mid_intervals', 'interv_mid_dur', 'tvals_mid',...
    'pvals_mid', 'dvals_mid')
    
%% 3.3.2 - Plot Figure 4b - effect size vs. interval centre

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/t_diff_mid_intervals')

figure, bar(interv_mid_dur,mean(dvals_mid))
hold on, errorbar(interv_mid_dur,mean(dvals_mid), ...
    std(dvals_mid)/sqrt(size(dvals_mid,1)))
anova_rm(dvals_mid)
[~,p,~,stat]=ttest(dvals_mid(:, interv_mid_dur==[-2])-dvals_mid(:, interv_mid_dur==[-1.5]));
fprintf('-1.5 > -2: t=%g, p=%g\n', stat.tstat, p)
[~,p,~,stat]=ttest(dvals_mid(:, interv_mid_dur==[-2.5])-dvals_mid(:, interv_mid_dur==[-1.5]));
fprintf('-1.5 > -2.5: t=%g, p=%g\n', stat.tstat, p)
[~,p,~,stat]=ttest(dvals_mid(:, interv_mid_dur==[-3])-dvals_mid(:, interv_mid_dur==[-1.5]));
fprintf('-1.5 > -3: t=%g, p=%g\n', stat.tstat, p)
    
%% 3.3.3 - Control analysis for effect size vs. interval centre
%%% do simulations of different centres of preseimtulus PEP interval
%%% show the same pattern as their impact on subsequent endorsements, which
%%% could imply that it is an outcome of different number of cardiac peaks

interv_mid_dur = -4.5:.5:-1.5;
nintervs = numel(interv_mid_dur);
interv_dur = 3;

effect_size = 5.5;
% based on the effect size difference between high- and
% low-arousal conditions in in Luft and Bhattacharya (2016)
% Figure 3 (in µV/m3)

nperms = 10000; % number of permutations
alpha = 0.05;

endorscutoff = .9;

probdetect = NaN(nintervs, 2);

for nb = 1:nintervs
    bef_aft_hep = [interv_mid_dur(nb)-interv_dur/2 ...
        interv_mid_dur(nb)+interv_dur/2];
    
    % first determine the mean number of PEP trials per subject
    % and per Valence (pleasant, unpleasant) x Endorse (Yes, No)
    % conditions. also determine the average noise in EEG
    % channels (defined as std of frequencies > 7 Hz) across
    % subjects, and the number of subjects for each threshold
    [hbtrials, noise_subj, val_subj] = hepsret_optim_inclus...
        (exp_hepsret, filename, ...
        bef_aft_hep, trange_hep, endorscutoff, pulsedelay);
    nsubj_valid = sum(val_subj);
    
    for pn = 1:2 % positive and negative
        % simulate the experiment with expected effect size and
        % average noise for each threshold
        probdetect(nb,pn) = hepsret_simulate_detect_prob...
            (effect_size, noise_subj, squeeze(round(hbtrials(:, pn, :))),...
            nsubj_valid, nperms, alpha);
    end
end

save('hepsret_datafiles/int_cent_power', 'nintervs', 'probdetect',...
    'interv_mid_dur', 'interv_dur', 'nperms', 'effect_size')

figure
plot(interv_mid_dur, probdetect(:,2), '-', 'LineWidth', 1.5)

%%% doesn't seem to be the case

%%
%% 
%%
%% 3.4.1 - Do post stimulus HEP analysis

bef_aft_hep = [0 5]; % look at post-stimulus HEPs

hepsret_hep(exp_hepsret, selchans, filename, ...
    smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr,...
    blsub_flag, endorscutoff,...
    do_plots);
    
%% 3.4.2 statistics  post stimulus HEPs for Pleasant-Self vs. Pleasant-Nonself
stat_type = 8; % for positive HEP
alpha = 0.025;

hepsret_stats(exp_hepsret, stat_type, neighbours, alpha);
%%% no significant clusters for pos hep for yes v. no

%% 3.4.3 statistics on post stimulus HEPs for UnPleasant-Self vs. UnPleasant-Nonself

stat_type = 9; % for negative HEP
hepsret_stats(exp_hepsret, stat_type, neighbours, alpha);
%%% no significant clusters for pos hep for yes v. no

%% 3.4.4 stats on post HEPs between pos and neg (not by endorsement)

stat_type = 10; % for plotting positive HEP

hepsret_stats(exp_hepsret, stat_type, neighbours, alpha);
