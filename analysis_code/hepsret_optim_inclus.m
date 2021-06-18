function [hbtrials, noise_subj, val_subj] = ...
    hepsret_optim_inclus(exp_hepsret, filename, ...
    bef_aft_hep, trange_hep, endorscutoff)
%
%
% to determine the optimal number of participants to include so that there
% are sufficient number of trials
%

if ~exist('bef_aft_hep', 'var')||isempty(bef_aft_hep)
    bef_aft_hep = [-.25 1.3];
end

chpktr = hepsret_changepeak;

% High pass filter
fs=250;
hp_filt_cutoff = 7;
filter_order = 3;
hpfilt_cutoff = hp_filt_cutoff/(fs/2);
[hpfilt_b, hpfilt_a] = butter(filter_order, hpfilt_cutoff, 'high');


nvalsubj = 0;

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/wordresp.mat')
hbchanname = 'BVP';

chanrem = hepsret_chanrem;

nsubj = exp_hepsret.nsubj;
val_subj = false(1, nsubj);

hb_endors = NaN(nsubj, 2, 2);
noise_subj = NaN(1, nsubj);

for ns = 1:nsubj
    
    subj_data = exp_hepsret.data(ns);
    
    subj_endors = all_word_resp(ns, :);
    %     if ns<=4
    %         subj_endors([26:30, 56:60]) = [];
    %         ntn = 25;
    %     else
    ntn = size(subj_endors,2)/2; % number of pleasant and unpleasant words
    %     end
    
    fprintf('\nsubj id: %s\n', subj_data.dir_name)
    pos = nansum(subj_endors(1:end/2));
    neg = nansum(1-subj_endors(end/2+1:end));
    fprintf('perc pos: %g, perc neg: %g\n', pos/ntn, neg/ntn)
    
    if subj_data.subj_valid_sret && ...
            subj_data.subj_code~=4169 && ... % 4169 has no physio data
            (nansum(subj_endors(1:end/2))/ntn<=endorscutoff &&...
            nansum(1-subj_endors(end/2+1:end))/ntn<=endorscutoff)
        
        cd(fullfile(exp_hepsret.session_dir, subj_data.dir_name))
        nvalsubj = nvalsubj+1;
        val_subj(ns) = true;

        posneg_wordinds = medriv_posneginds(subj_data.subj_valid_sret);
        
        EEG = pop_loadset(['Data/eeglab/' filename '.set']);
        
        %%% select only the eeg channels
        chansel = false(1, EEG.nbchan);
        chanlabs = {EEG.chanlocs.labels};
        keepchans = ~get_channels_from_labels(chanlabs, chanrem);
        keepchans(32:end) = false;
        
        hbchan = find(strcmp(chanlabs, hbchanname));
        
        if ~isempty(hbchan)
            
            %%% get hep
            [~,hb_posneg] = sret_hep_analysis(EEG, chansel, hbchan, ...
                posneg_wordinds, bef_aft_hep, trange_hep,...
                false, false, false,...
                chpktr([chpktr.subj]==subj_data.subj_code));
            % epochs - num of words x before after 0 x channels x time
            
            subj_endors = reshape(subj_endors, [ntn 2])';
            for npp = 1:2 % heps pre post stimuulus
                for pn = 1:2 % pos neg
                 
                    %%% count the number of endorsements
                    hb_endors(nvalsubj, ...
                        pn, 1) = sum(hb_posneg(subj_endors(pn, :)==1, pn));
                    hb_endors(nvalsubj, ...
                        pn, 2) = sum(hb_posneg(subj_endors(pn, :)==0, pn));
                end
            end
            
            % high pass filter the data above theta frequencies
            eegdata = double(EEG.data(keepchans,:));
            nchans = size(eegdata,1);
            for nch = 1:nchans
                eegdata(nch,:) = filtfilt(hpfilt_b, hpfilt_a, eegdata(nch,:));
            end

            noise_subj(ns) = mean(std(eegdata, 0, 2));

        end
        
    end
end

hbtrials = hb_endors(~isnan(hb_endors(:,1,1)),:,:);
% plot the mean and stderr of number of trials for each condition
%  = squeeze(mean(hb_endors));
% shb = squeeze(std(hb_endors)/sqrt(nvalsubj));
% figure, errorbar(meantrials, shb)
noise_subj = noise_subj(~isnan(noise_subj));
% fprintf('\naverage eeg noise = %g\n', mean(noise_subj))

end