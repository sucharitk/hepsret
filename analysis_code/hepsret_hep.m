function hepsret_hep(exp_hepsret, selchans, filename, ...
    smooth_erp, bef_aft_hep, trange_hep, bef_aft_erp, bef_aft_hr, ...
    blsub_flag, endorscutoff, shuffle_peaks, do_plots)


if ~exist('bef_aft_hep', 'var')||isempty(bef_aft_hep)
    bef_aft_hep = [-.25 1.3];
end
if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end
if ~exist('shuffle_peaks', 'var')||isempty(shuffle_peaks), shuffle_peaks = false; end
if ~exist('do_plots', 'var')||isempty(do_plots), do_plots = false; end

nselchan = numel(selchans);
nvalsubj = 0;
eeg_mode = 1;

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/wordresp.mat')

hbchanname = exp_hepsret.hbchanname;
hrchanname= exp_hepsret.hrchanname;

% channels to exclude from analysis
chanrem = hepsret_chanrem;

srate = 250;
nhepinds = ceil(diff(trange_hep*srate)+.5);
nerpinds = round(diff(bef_aft_erp*srate))+1;
nhrinds =  round(diff(bef_aft_hr*srate))+1;

nsubj = exp_hepsret.nsubj;

% manually add or remove a few incorrectly detected peaks
chpktr = hepsret_changepeak;

sret_hep = NaN(12, 2, 2, 32, nhepinds); % do not split HEP by endorsement
sret_hep_endors = NaN(12, 2, 2, 2, 32, nhepinds); % split HEP by endorsement
sret_erp_endors = NaN(12, 2, 2, 32, nerpinds);  % erp 
sret_hr_endors = NaN(12, 2, 2, nhrinds);  % hr-erp
total_hb = 0;

for ns = 1:nsubj
    
    subj_data = exp_hepsret.data(ns);
    
    subj_endors = all_word_resp(ns, :);

    ntn = size(subj_endors,2)/2; % number of pleasant and unpleasant words

    
    fprintf('\nsubj id: %s\n', subj_data.dir_name)
    pos = nansum(subj_endors(1:end/2));
    neg = nansum(1-subj_endors(end/2+1:end));
    fprintf('perc pos: %g, perc neg: %g\n', pos/ntn, neg/ntn)
    
    if subj_data.subj_valid_sret &&... % 4163's heart data is basically noise
            subj_data.subj_code~=4169 &&... %4163's  has no heart data 
            (pos/ntn<=endorscutoff &&...
            neg/ntn<=endorscutoff)
        
        cd(fullfile(exp_hepsret.session_dir, subj_data.dir_name))
        nvalsubj = nvalsubj+1;
        
        posneg_wordinds = medriv_posneginds(subj_data.subj_valid_sret);
        
        EEG = pop_loadset(['Data/eeglab/' filename '.set']);
        
        %%% select the channel indices
        chanlabs = {EEG.chanlocs.labels};
        if isempty(selchans)
            if ~exist('chanlocs', 'var')
                chanlocs = EEG.chanlocs;
                chanlocs = chanlocs(1:32);
                nselchan = numel(chanlocs);
                selchans = chanlabs;
            end
        end
        eegchans = 1:nselchan;
        [chansel, chansavail] = rearrange_channels(chanlabs, selchans,...
            nselchan);
        
        hbchan = find(strcmp(chanlabs, hbchanname));
        hrchan = find(strcmp(chanlabs, hrchanname));
        
        if ~isempty(hbchan)
            
            hdata = EEG.data(hbchan,:);
            hr = hr_calc(hdata', srate);
            ntrials(nvalsubj, :) = hr*diff(bef_aft_hep)*[ntn-neg ntn-pos]/60;
            fprintf('ntrials: pos=%g, neg=%g\n', ntrials(nvalsubj, 1), ...
                ntrials(nvalsubj, 2))
            %%% get hep
            [epochs,hb_posneg] = sret_hep_analysis(EEG, chansel, hbchan, ...
                posneg_wordinds, bef_aft_hep, trange_hep,...
                smooth_erp, blsub_flag, shuffle_peaks, ...
                chpktr([chpktr.subj]==subj_data.subj_code), do_plots);
            % epochs - num of words x before after 0 x channels x time
            
            total_hb = total_hb + sum(hb_posneg(:));
            
            %%% get erp
            epochs2 = sret_erp_analysis(EEG, eeg_mode, posneg_wordinds, ...
                bef_aft_erp, eegchans, smooth_erp, blsub_flag);
            
            %%% get Heartrate
            % first smooth the heart rate
            EEG.data(hrchan, :) = smooth(EEG.data(hrchan, :), 100); % smooth the heartrate
            EEG.data(hrchan, :) = EEG.data(hrchan, :) - mean(EEG.data(hrchan, :)); % demean the heartrate data
            epochs3 = sret_erp_analysis(EEG, eeg_mode, posneg_wordinds, ...
                bef_aft_hr, hrchan, smooth_erp, false); % do again without baseline subtraction for the heartrate data

            if do_plots
                figure(2)
                clf, hold on
            end
            
            subj_endors = reshape(subj_endors, [ntn 2])';
            for npp = 1:2 % heps pre post stimuulus
                for pn = 1:2 % pos neg
                    sret_hep(nvalsubj, ...
                        pn, npp, chansavail, :) = ...
                        squeeze(nanmean(epochs{pn}(:, npp, :, :), 1));
                    if do_plots, plot( squeeze(nanmean(nanmean(epochs{pn}...
                            (:, npp, :, :), 1),3))'), end
                    
                    %%% split the heps by endorsements
                    sret_hep_endors(nvalsubj, ...
                        pn, npp, ...
                        1, chansavail, :) = ...
                        squeeze(nanmean(epochs{pn}(subj_endors(pn, :)==1, npp, ...
                        :, :), 1)); % yes
                    sret_hep_endors(nvalsubj, ...
                        pn, npp, ...
                        2, chansavail, :) = ...
                        squeeze(nanmean(epochs{pn}(subj_endors(pn, :)==0, npp, ...
                        :, :), 1)); % no
                    
                end
            end
            
            for pn = 1:2 % pos neg
                %%% split the erps by endorsements
                sret_erp_endors(nvalsubj, pn, 1, eegchans, :) = permute(nanmean(...
                    epochs2{pn}(subj_endors(pn, :)==1, eegchans, :),1), [2:3 1]); % yes endorsements
                sret_erp_endors(nvalsubj, pn, 2, eegchans, :) = permute(nanmean(...
                    epochs2{pn}(subj_endors(pn, :)==0, eegchans, :),1), [2:3 1]); % no endorsements
                
                %%% split hr erp by endorsements
                sret_hr_endors(nvalsubj, pn, 1, :) = permute(nanmean(...
                    epochs3{pn}(subj_endors(pn, :)==1, hrchan, :),1), [2:3 1]); % yes endorsements
                sret_hr_endors(nvalsubj, pn, 2, :) = permute(nanmean(...
                    epochs3{pn}(subj_endors(pn, :)==0, hrchan, :),1), [2:3 1]); % no endorsements
            end
        end
        
    end
end

calc_num_manual_peaks = false;
if calc_num_manual_peaks
    % if true, then calculate the percentage of hb peaks that were added or
    % removed manually 
    nch_pk = numel(chpktr);
    pk_add_sub = [0 0];
    for ncp = 1:nch_pk
        pp = abs(chpktr(ncp).cp(:,1));
        pk_add_sub(1) = pk_add_sub(1) + sum(pp<100);
        pk_add_sub(2) = pk_add_sub(2) + sum(pp>100);
    end
    fprintf('perc peaks manually added=%g, removed=%g', pk_add_sub(1)/total_hb,...
        pk_add_sub(2)/total_hb)
    
end
sret_hep_endors(sret_hep_endors==0) = NaN;

labchans = {chanlocs.labels};
keepchans = ~get_channels_from_labels(labchans, chanrem);
labchans = labchans(keepchans);

srethep = sret_hep_endors(:, :, :, :, keepchans, :);

hrep = squeeze(srethep(:,:,1,:,:,:)); %  save the prestimulus value
hrep_post = squeeze(srethep(:,:,2,:,:,:)); %  save the poststimulus value

hrep_post_all = squeeze(sret_hep(:,:,2,keepchans,:)); % save the poststimulus hep unseparated by endorsement

erp_sret = sret_erp_endors(:,:,:,keepchans,:); % exclude the peripheral channels
hr_sret = sret_hr_endors(:, :, :, :);

cd(exp_hepsret.session_dir)

if shuffle_peaks
    save('hepsret_datafiles/hrep_shuffle','hrep',...
        'labchans','srate', 'trange_hep')
else
    save('hepsret_datafiles/hrep','hrep','hrep_post', 'hrep_post_all',...
        'labchans','srate', 'trange_hep',...
        'erp_sret', 'hr_sret', 'bef_aft_erp', 'bef_aft_hr')
end

end