function [pc1,pc2] = hepsret_decode_selfendors(exp_hepsret, stats, filename,...
    bef_aft_hep, trange_hep, bef_aft_hr,...
    blsub_flag, endorscutoff, nshuf, learn_method)
% 
% hepsret_decode_selfendors(exp_hepsret, stats, filename,...
%     bef_aft_hep, trange_hep, bef_aft_hr,...
%     blsub_flag, endorscutoff, nshuf, learn_method)
%
% trial-by-trial decoding of self-endorsements from pre-stimulus HEPs using
% leave-1-out logistic regression / svm
%


if ~exist('bef_aft_hep', 'var')||isempty(bef_aft_hep)
    bef_aft_hep = [-.25 1.3];
end
if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end

eeg_mode = 1;
nvalsubj = 0;

cd(exp_hepsret.session_dir)
load('hepsret_datafiles/wordresp.mat')

hbchanname = 'BVP';

chanrem = hepsret_chanrem;
statclust = stats.posclusterslabelmat==1;
statclust = stats.posclusterslabelmat<5 & stats.posclusterslabelmat>0;

srate = 250;
nhepinds = ceil(diff(trange_hep*srate)+.5);
nhrinds =  round(diff(bef_aft_hr*srate))+1;
thep = linspace(trange_hep(1), trange_hep(2), nhepinds);
thr = linspace(bef_aft_hr(1), bef_aft_hr(2), nhrinds);

nsubj = exp_hepsret.nsubj;

catvar = true;

chpktr = hepsret_changepeak;
pc1 = NaN(12, nshuf);
pc2 = NaN(12, nshuf);

for ns = 1:nsubj
    
    subj_data = exp_hepsret.data(ns);
    
    subj_endors = all_word_resp(ns, :);
    ntn = size(subj_endors,2)/2; % number of pleasant and unpleasant words
    
    if subj_data.subj_valid_sret && subj_data.subj_code~=4163 &&... % 4163's heart data is basically noise
            subj_data.subj_code~=4169 &&... %4163's  has no heart data
            (nansum(subj_endors(1:end/2))/ntn<=endorscutoff &&...
            nansum(1-subj_endors(end/2+1:end))/ntn<=endorscutoff)
        
        cd(fullfile(exp_hepsret.session_dir, subj_data.dir_name))
        nvalsubj = nvalsubj+1;
        
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
        eegchans = 1:nselchan;
        chansel = rearrange_channels(chanlabs, selchans,...
            nselchan);
        labchans = {chanlocs.labels};
        keepchans = ~get_channels_from_labels(labchans, chanrem);
        
        
        hbchan = find(strcmp(chanlabs, hbchanname));
        %         hrchan = find(strcmp(chanlabs, hrchanname));
        
        if ~isempty(hbchan)
            
%             hdata = EEG.data(hbchan,:);
            %             hr = hr_calc(hdata', srate);
            
            %%% get hep
            [~,~, indiv_hep, indiv_trialnums] = sret_hep_analysis...
                (EEG, chansel, hbchan, ...
                posneg_wordinds, bef_aft_hep, trange_hep,...
                [], blsub_flag, [], ...
                chpktr([chpktr.subj]==subj_data.subj_code));
            % epochs - num of words x before after 0 x channels x time
            
            hepd = cell(1, 2);
            subj_endors = reshape(subj_endors, [ntn 2])';
            

            for pn = 1:2 % pos neg
                nthep = numel(indiv_trialnums{pn});
                for nt = 1:nthep
                    h2 = squeeze(indiv_hep{pn}(nt, keepchans, thep>=0));
                    
                    %%% for using only a selected number of time points for
                    %%% the significant channels
                    hepd{pn}(nt, :) = h2(statclust);
                    
                    %%% for using all time points of the significant
                    %%% channels
                    %                     h2 = h2(logical(sum(statclust,2)),:);
                    %                     hepd{pn}(nt, :) = h2(:);
                end
            end
            
            % expand the endorsements to match the number of heartbeat
            % trials
            trialvalsout = expand_trial2hep(subj_endors(pn, :),...
                indiv_trialnums{pn});
            
            for pn = 2:2 % pos neg
                % do logistic regression here
                hepd_pn = hepd{pn};
                
                e1 = find(trialvalsout==1); % yes endorsements
                e0 = find(trialvalsout==0); % no endorsements
                n1 = numel(e1); n0 = numel(e0);
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
                
                num_feat = size(hepd_pn,2);
                
                X = NaN(tot_test_trials, num_feat, nshuf);
                Y = NaN(tot_test_trials, nshuf);
                X2 = X;
                
                if catvar
                    % if to do analyses by making the output variable
                    % categorical
                    trialvalsout = categorical(trialvalsout);
                     Y = categorical(Y);
                end
                
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
                    pc1(nvalsubj, nsh)=sum(Y(:,nsh)==fit)/tot_test_trials;
                    
                    % fit the model with the shuffled data
                    fit = leave1out_logreg(X2(:,:,nsh),Y(:,nsh), learn_method);
                    %                     pc2(nvalsubj, nsh) = sum(1-abs(Y-fit))/tot_test_trials;
                    pc2(nvalsubj, nsh) = sum(Y(:,nsh)==fit)/tot_test_trials;
                    
                end
              
                
                [nanmean(pc1,2) nanmean(pc2,2)]
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
%
end

function fit = leave1out_logreg(X,Y, method)

if ~exist('method', 'var'), method =1; end

ntr = size(X,1);
indices = 1:ntr;
fit = NaN(ntr,1);
if iscategorical(Y)
    fit = categorical(fit);
end



switch method
    case {1, 2, 3}
        for nt = 1:ntr
            test = indices==nt; train = ~test;
            
            switch method
                case 1
                    
                    b = glmfit(X(train,:),Y(train),...
                        'binomial','logit'); % Logistic regression
                    fit(nt) = round(glmval(b,X(test,:),'logit')');
                    
                case 2
                    
                    Mdl = fitclinear(X(train,:),Y(train), 'Learner', 'logistic');
                    fit(nt) = predict(Mdl,X(test,:));
                    
                    
                case 3
                    
                    Mdl = fitclinear(X(train,:),Y(train), 'Learner', 'svm');
                    fit(nt) = predict(Mdl,X(test,:));
                    
                    
                case 4
                    % better do 1, 4 takes longer
                    
                    b = fitglm(X(train,:),Y(train), 'link', 'logit',...
                        'Distribution', 'binomial');
                    fit(nt) = round(predict(b, X(test,:)));
                    
            end
        end
    case 4
end
end
