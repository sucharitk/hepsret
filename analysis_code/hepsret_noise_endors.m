function noisehep_subj = hepsret_noise_endors(exp_hepsret, ...
    stats, filename,...
    bef_aft_hep, trange_hep,...
    blsub_flag, endorscutoff)

if ~exist('bef_aft_hep', 'var')||isempty(bef_aft_hep)
    bef_aft_hep = [-.25 1.3];
end
if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end

nvalsubj = 0;

cd(exp_hepsret.session_dir)
load('all_word_resp.mat')

hbchanname = 'BVP';

chanrem = hepsret_chanrem;
statclust = stats.posclusterslabelmat==1;

srate = 250;
nhepinds = ceil(diff(trange_hep*srate)+.5);
thep = linspace(trange_hep(1), trange_hep(2), nhepinds);

nsubj = exp_hepsret.nsubj;

chpktr = hepsret_changepeak;

noisehep_subj = NaN(nsubj, 2, 2); % subj x (pos neg both) x (mean var)
bias_subj = NaN(nsubj, 2); % subj x (pos neg)
val_subj = false(1, nsubj);

for ns = 1:nsubj
    
    subj_data = exp_hepsret.data(ns);
    
    subj_endors = all_word_resp(ns, :);
    if ns<=4
        subj_endors([26:30, 56:60]) = [];
        ntn = 25;
    else
        ntn = 30;
    end
    
    bias_subj(ns, 1) = nansum(subj_endors(1:end/2))/ntn;
    bias_subj(ns, 2) = nansum(1-subj_endors(end/2+1:end))/ntn;
    
    if subj_data.subj_valid_sret && subj_data.subj_code~=4163 &&... % 4163's heart data is basically noise
            subj_data.subj_code~=4169 &&... %4163's  has no heart data
            (bias_subj(ns, 1)<=endorscutoff &&...
            bias_subj(ns, 2)<=endorscutoff)
        
        cd(fullfile(exp_hepsret.session_dir, subj_data.dir_name))
        nvalsubj = nvalsubj+1;
        val_subj(ns) = true;
        
        posneg_wordinds = medriv_posneginds(subj_data.subj_valid_sret);
        
        EEG = pop_loadset(['Data/eeglab/' filename '.set']);
        
        %%% select the channel indices
        chanlabs = {EEG.chanlocs.labels};
        if ~exist('selchans', 'var')||isempty(selchans)
            if ~exist('chanlocs', 'var')
                chanlocs = EEG.chanlocs;
                chanlocs = chanlocs(1:32);
                nselchan = numel(chanlocs);
                selchans = chanlabs;
            end
        end
        
        [chansel, ~] = rearrange_channels(chanlabs, selchans,...
            nselchan);
        labchans = {chanlocs.labels};
        keepchans = ~get_channels_from_labels(labchans, chanrem);
        
        
        hbchan = find(strcmp(chanlabs, hbchanname));
        
        if ~isempty(hbchan)
            
            %%% get hep
            [~,~, indiv_hep, indiv_trialnums] = sret_hep_analysis...
                (EEG, chansel, hbchan, ...
                posneg_wordinds, bef_aft_hep, trange_hep,...
                [], blsub_flag, [], ...
                chpktr([chpktr.subj]==subj_data.subj_code));
                        
            hepd = cell(1, 2);
            
            subj_endors = reshape(subj_endors, [ntn 2])';

            mh = hepd;
            for pn = 1:2 % pos neg
                nthep = numel(indiv_trialnums{pn});
                for nt = 1:nthep
                    h2 = squeeze(indiv_hep{pn}(nt, keepchans, thep>=0));
                    hepd{pn}(nt, :) = h2(statclust);
                    
                end
                mh{pn} = mean(hepd{pn},2);
                trialvalsout = expand_trial2hep(subj_endors(pn, :),...
                    indiv_trialnums{pn});
                
                noisehep_subj(ns, pn, 1) = std(mh{pn}(trialvalsout==1)); % yes
                noisehep_subj(ns, pn, 2) = std(mh{pn}(trialvalsout==0)); % no

            end
            
            
        end
        
    end
end

noisehep_subj = noisehep_subj(val_subj,:,:);
end


function trialvalsout = expand_trial2hep(trialvalsin, trialnumsto)
%
% takes individual trial and then expands them to the number of trials
% based on the HEPs
%

ntrialto = numel(trialnumsto);
trialvalsout = NaN(1, ntrialto);
for ntt = 1:ntrialto
    trialvalsout(ntt) = trialvalsin(trialnumsto(ntt));
end
%
end

