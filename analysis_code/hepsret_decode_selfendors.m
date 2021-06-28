function [pc1,pc2] = hepsret_decode_selfendors(exp_hepsret, stats, filename,...
    bef_aft_hep, trange_hep, bef_aft_hr,...
    blsub_flag, endorscutoff, shuf_chans, nshuf, learn_method)
%
% hepsret_decode_selfendors(exp_hepsret, stats, filename,...
%     bef_aft_hep, trange_hep, bef_aft_hr,...
%     blsub_flag, endorscutoff, nshuf, learn_method)
%
% trial-by-trial decoding of self-endorsements from pre-stimulus HEPs using
% leave-1-out logistic regression / svm
%

subselect = false;

if ~exist('bef_aft_hep', 'var')||isempty(bef_aft_hep)
    bef_aft_hep = [-.25 1.3];
end
if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end

nvalsubj = 0;

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/wordresp.mat')
load('chanlocs.mat')
chanlab32 = {chanlocs32.labels};
neegchans = numel(chanlab32);

hbchanname = 'BVP';

plotchan = {'Fz'};

chanrem = hepsret_chanrem;
sigclust = find([stats.posclusters.prob]<.05);
% statclust = stats.posclusterslabelmat==1;
statclust = stats.posclusterslabelmat<=sigclust(end) & stats.posclusterslabelmat>0;

srate = 250;
nhepinds = ceil(diff(trange_hep*srate)+.5);
% nhrinds =  round(diff(bef_aft_hr*srate))+1;
thep = linspace(trange_hep(1), trange_hep(2), nhepinds);
% thr = linspace(bef_aft_hr(1), bef_aft_hr(2), nhrinds);

nsubj = exp_hepsret.nsubj;


chpktr = hepsret_changepeak;

pc1 = NaN(12,1);
if nargout>1
    % do the shuffled analysis
    pc2 = NaN(12, nshuf);
    permut_pc = true; % do permutation analyses
    
else
    % if only one output is expected, don't do the shuffling, only return
    % the pc correct
    permut_pc = false;
    plotchan = {'Oz'};
    chans_to_shuf = {'Oz', 'O1', 'O2'};
    
end


fig_handle=figure;subplot(531)

for ns = 1:nsubj
    
    subj_data = exp_hepsret.data(ns);
    
    subj_endors = all_word_resp(ns, :);
    ntn = size(subj_endors,2)/2; % number of pleasant and unpleasant words
    
    if subj_data.subj_valid_sret && ...
            subj_data.subj_code~=4169 &&... %4163's  has no heart data
            (nansum(subj_endors(1:end/2))/ntn<=endorscutoff &&...
            nansum(1-subj_endors(end/2+1:end))/ntn<=endorscutoff)
        
        cd(fullfile(exp_hepsret.session_dir, subj_data.dir_name))
        nvalsubj = nvalsubj+1;
        
        fprintf('ns=%g, nval=%g - %g\n', ns, nvalsubj, subj_data.subj_code)
        
        posneg_wordinds = medriv_posneginds(subj_data.subj_valid_sret);
        
        EEG = pop_loadset(['Data/eeglab/' filename '.set']);
        
        %%% select the channel indices
        chanlabs_curds = {EEG.chanlocs.labels};
        
        
        eegchans = rearrange_channels(chanlabs_curds, chanlab32,...
            neegchans);
        
        keepchans = ~get_channels_from_labels(chanlabs_curds(eegchans),...
            chanrem); % remove the excluded eeg channels from the current dataset
        
        
        hbchan = find(strcmp(chanlabs_curds, hbchanname));
        
        if ~isempty(hbchan)
            
            %%% get hep
            [~,~, indiv_hep, indiv_trialnums] = sret_hep_analysis...
                (EEG, eegchans, hbchan, ...
                posneg_wordinds, bef_aft_hep, trange_hep,...
                [], blsub_flag, [], ...
                chpktr([chpktr.subj]==subj_data.subj_code));
            % epochs - num of words x before after 0 x channels x time
            
            hepd = cell(1, 2);
            subj_endors = reshape(subj_endors, [ntn 2])';
            
            for pn = 2:2 % pos neg
                
                nthep = numel(indiv_trialnums{pn});
                
                for nt = 1:nthep
                    %                     h2 = squeeze(indiv_hep{pn}(nt, keepchans, thep>=0));
                    h2 = squeeze(indiv_hep{pn}(nt, keepchans, thep>=0));
                    
                    %%% for using only a selected number of time points for
                    %%% the significant channels
                    
                    if shuf_chans
                        % switch place of Fz, FC1 and FC2 with occipital
                        % channels Oz, O1, O2
                        shufchaninds = ...
                            get_channels_from_labels(...
                            chanlabs_curds(keepchans),chans_to_shuf);
                        
                        chans_to_swapout = logical(mean(statclust,2))';
                        
                        h2(chans_to_swapout,:) = h2(shufchaninds,:);
                        
                    end
                    
                    h2 = (h2(statclust));
                    
                    %                     sigchans = find(mean(statclust,2));
                    %                     nsigchan = numel(sigchans);
                    %                     h3 = [];
                    %                     for nsc = 1:nsigchan
                    %                         h3 = [h3 ...
                    %                             mean(h2(sigchans(nsc), statclust(sigchans(nsc),:)))];
                    %                     end
                    %                     h2=h3;
                    
                    hepd{pn}(nt, :) = h2(:);
                end
            end
            
            
            % expand the endorsements to match the number of heartbeat
            % trials
            fzc = get_channels_from_labels(chanlabs_curds, plotchan);
            trialvalsout = expand_trial2hep(subj_endors(pn, :),...
                indiv_trialnums{pn});
            
            figure(fig_handle)
            subplot(5,3,nvalsubj),hold off,
            plot(permute(mean(indiv_hep{pn}(trialvalsout==1,fzc,:),1), [3 2 1]))
            hold on
            plot(permute(mean(indiv_hep{pn}(trialvalsout==0,fzc,:),1), [3 2 1]))
            
            for pn = 2:2 % pos neg
                
                % do logistic regression here
                hepd_pn = hepd{pn};
                
                if size(h2,1)>1
                    chanbytime = repmat(1:size(statclust,1),[size(statclust,2) 1])';
                    chanbytime = chanbytime(statclust);
                    
                    hepd_pn = standnorm(hepd_pn, 1, chanbytime);
                    
                else
                    % standardise/normalise the data
                    hepd_pn = standnorm(hepd_pn, 1);
                    
                end
                
                e1 = find(trialvalsout==1); % yes endorsements
                e0 = find(trialvalsout==0); % no endorsements
                n1 = numel(e1); n0 = numel(e0);
                num_feat = size(hepd_pn,2);
                
                
                if  subselect
                    if n0 > n1
                        % num trials per condition is for the ones with the
                        % fewer trials
                        ntp_cond = n1;
                        
                        % num of shuffles of the condition with more trials
                        nshuf_cond = n0;
                        
                        % endorsements that need to be subselected
                        esel = e0;
                        
                        % endorsements that need to be all selected
                        eall = e1;
                    else
                        
                        ntp_cond = n0;
                        nshuf_cond = n1;
                        esel = e1;
                        eall = e0;
                        
                    end
                    % total test trials
                    tot_test_trials = ntp_cond*2;
                    
                    
                    X = NaN(tot_test_trials, num_feat, nshuf);
                    Y = NaN(tot_test_trials, nshuf);
                    X2 = X;
                    
                    %                     if catvar
                    %                         % if to do analyses by making the output variable
                    %                         % categorical
                    %                         trialvalsout = categorical(trialvalsout);
                    %                         Y = categorical(Y);
                    %                     end
                    
                    for nsh = 1:nshuf
                        
                        % shuffle the trials that need to subselected
                        esel_shuf_all = esel(randperm(nshuf_cond));
                        % subselect only the number of trials in the condition
                        % with fewer number of trails
                        esel_shuf = esel_shuf_all(1:ntp_cond);
                        % append it to the condition with the fewere number of
                        % trials
                        trialsel = [eall esel_shuf];
                        
                        % set X and Y
                        X(:,:,nsh) = hepd_pn(trialsel,:);
                        Y(:,nsh) = trialvalsout(trialsel);
                        
                        
                        % shuffle X to get a measure of chance percent correct
                        shuf_control = randperm(ntp_cond*2);
                        X2(:,:,nsh) = X(shuf_control, :,nsh);
                        
                    end
                    
                    parfor nsh = 1:nshuf
                        
                        % first fit the model with the actual data
                        fit = leave1out_logreg(X(:,:,nsh),Y(:,nsh), learn_method);
                        %                     pc1(nvalsubj, nsh)=sum(1-abs(Y-fit))/tot_test_trials;
                        % calculate and save percent correct
                        pc1(nvalsubj, nsh) = sum(Y(:,nsh)==fit)/tot_test_trials;
                        
                        % fit the model with the shuffled data
                        fit = leave1out_logreg(X2(:,:,nsh),Y(:,nsh), learn_method);
                        %                     pc2(nvalsubj, nsh) = sum(1-abs(Y-fit))/tot_test_trials;
                        pc2(nvalsubj, nsh) = sum(Y(:,nsh)==fit)/tot_test_trials;
                        
                    end
                    
                    figure(fig_handle)
                    subplot(5,3,nvalsubj),
                    title(sprintf('ns=%g, pc1=%g, pc2=%g', ns,...
                        nanmean(pc1(nvalsubj,:),2), nanmean(pc2(nvalsubj,:),2)))
                    
                    
                else
                    
                    
                    X = hepd_pn;
                    Y = trialvalsout';
                    nany = isnan(Y);
                    X(nany,:) =[];
                    Y(nany) = [];
                    tot_test_trials = n1+n0;
                    
                    fit = leave1out_logreg(X,Y, learn_method);
                    
                    pc1(nvalsubj) = sum(Y==fit)/tot_test_trials;
                    
                    title(sprintf('ns=%g, pc1=%g', ns,...
                        pc1(nvalsubj)))
                    
                    
                    if permut_pc
                        % don't do this permutation part if there is no
                        % second return argument - this will happen if we
                        % are shuffling the channels to compare the
                        % effect with a control channel cluster
                        
                        X2 = NaN(tot_test_trials, num_feat, nshuf);
                        for nsh = 1:nshuf
                            % shuffle X to get a measure of chance percent correct
                            X2(:,:,nsh) = X(randperm(tot_test_trials), :);
                            
                        end
                        
                        
                        parfor nsh = 1:nshuf
                            % fit the model with the shuffled data
                            fit2 = leave1out_logreg(X2(:,:,nsh),Y, learn_method);
                            %                     pc2(nvalsubj, nsh) = sum(1-abs(Y-fit))/tot_test_trials;
                            pc2(nvalsubj, nsh) = sum(Y==fit2)/tot_test_trials;
                        end
                        
                        figure(fig_handle)
                        subplot(5,3,nvalsubj),
                        title(sprintf('ns=%g, pc1=%g, pc2=%g', ns,...
                            pc1(nvalsubj), nanmean(pc2(nvalsubj,:),2)))
                        
                    end
                    
                end
            end
            
            
        end
        
    end
end


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

end

function sndata = standnorm(data, stand_or_norm, featids)

if size(data,2)==1
    switch stand_or_norm
        case 1
            % standardise the data
            sndata = (data-mean(data))/std(data);
        case 2
            % normalise the data
            sndata = (data-max(data))/(max(data)-min(data));
    end
else
    numfeat = unique(featids);
    nfeats = numel(numfeat);
    sndata = data;
    for nf = 1:nfeats
        selfeat = featids==numfeat(nf);
        seldata = data(:, selfeat);
        switch stand_or_norm
            case 1
                % standardise the data                
                sndata(:,selfeat) = ...
                    (data(:,selfeat)-mean(seldata(:)))/std(seldata(:));
            case 2
                % normalise the data
                sndata(:,selfeat) = ...
                    (data(:,selfeat)-max(seldata(:)))/...
                    (max(seldata)-min(seldata));
        end
        
    end
end
end
