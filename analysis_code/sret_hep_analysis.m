function [epochs,hb_posneg, indiv_hep, indiv_trialnums] = sret_hep_analysis...
    (EEG, selchans, hchan, ...
    posneg_wordinds, bef_aft, trange_hep,...
    smooth_erp, blsub_flag, shuffle_peaks, chpk, do_plots)

%
%
%
% indiv_hep: returns the individual heps not averaged by trial, trialnum:
% returns the individual trial numbers corresponding to the individual HEPs

if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end
if ~exist('smooth_erp', 'var')||isempty(smooth_erp), smooth_erp = 0; end
if ~exist('pulsedelay', 'var')||isempty(pulsedelay), pulsedelay = 0; end
if ~exist('shuffle_peaks', 'var')||isempty(shuffle_peaks), shuffle_peaks = false; end

if nargout>2, save_indiv_hep = true; else, save_indiv_hep = false; end

trigsvals = cell(1, 2);
artif_trigs = {'rej1', 'rej2'};

for pn = 1:2
    tinds = posneg_wordinds{pn};
    ntind = numel(tinds);
    winds = cell(1, ntind);
    for nt = 1:ntind
        winds{nt} = add_actichamp_S(tinds(nt));
    end
    trigsvals{pn} = winds;
end
srate = EEG.srate;
indrange_hep = trange_hep*srate;
pulsedelay = pulsedelay*srate;
pkrange = [.1 .2]*EEG.srate; % duration after the detected peak to find the actual peak

neegchans = numel(selchans);
%%% get epochs corresponding to pos and neg words


baseinds = 1:(-trange_hep(1)*srate);

nhepinds = ceil(diff(indrange_hep)+0.5);

epochs = cell(1, 2);
if save_indiv_hep
    indiv_hep = epochs; indiv_trialnums = epochs;
end
% val_trials = cell(1, 2);
if do_plots
    figure(1)
    nerpinds = ceil(diff(bef_aft)*srate+0.5);
    tpl = linspace(bef_aft(1),bef_aft(2),nerpinds);
end
for pn = 1:2
    tinds = trigsvals{pn};
    ntind = numel(tinds);
    
    
    ep1 = NaN(ntind, 2, neegchans, nhepinds); % second dim is for pre and post 0
    ep3 = []; trialnums = [];
    
    for nt = 1:ntind
        
        trialnum = str2double(tinds{nt}(end-1:end));
        
        %         if isempty(noistrial)||~any(trialnum==[noistrial.tr])
        try
            % when using windows > 5 sec this will throw an error for the
            % very first trial because the data only has 5 seconds
            % preceding the stimulus for the first trial
            ep = Get_Epochs(EEG, tinds{nt}, [], ...
                floor(bef_aft*srate), [], artif_trigs);
            nerpinds = round(diff(bef_aft*srate))+1;
            tt = linspace(bef_aft(1), bef_aft(2), nerpinds);
            if bef_aft(2)>=0
                zind = find(tt==0);
            else
                zind = find(tt==bef_aft(2));
            end
            
        catch
            % in case of this error use a smaller (5 sec window for the
            % first trial
            tempbefaft = [-5 bef_aft(2)];
            ep = Get_Epochs(EEG, tinds{nt}, [], ...
                floor(tempbefaft*srate), [], artif_trigs);
            nerpinds = round(diff(tempbefaft*srate))+1;
            tt = linspace(tempbefaft(1), tempbefaft(2), nerpinds);
            if bef_aft(2)>=0
                zind = find(tt==0);
            else
                zind = find(tt==bef_aft(2));
            end
            
        end
        
        bhep = []; ahep = []; hb = zeros(1, 2);
        
        if do_plots, clf, hold on, title(sprintf('pn=%g, tr=%g',pn,nt)), end
        if ~isempty(ep)
            ep2 = ep{1};
            hev = ep2(hchan,:);
            if do_plots==1, plot(tpl,hev), end
            
            if ~isempty(hev)
                rrp = rr_peaks(hev', srate);
                
                if ~isempty(chpk)
                    
                    % for the few trials where some peaks need to be
                    % removed or added
                    ncp = numel(chpk);
                    for nchp = 1:ncp
                        if chpk(nchp).pn==pn && chpk(nchp).tn==nt
                            nchng = size( chpk(nchp).cp,1);
                            for ncg = 1:nchng
                                chtr = chpk(nchp).cp(ncg, :);
                                if abs(chtr(1))>=100
                                    
                                    % removing peaks are indicated by
                                    % adding 100 to the peak time range
                                    rrp(rrp>=chtr(1)+100-bef_aft(1) & ...
                                        rrp<=chtr(2)+100-bef_aft(1))...
                                        = [];
                                else
                                    
                                    % addpeak
                                    chtr = (chtr-tt(1))*srate;
                                    if all(chtr>0) && all(chtr<nerpinds)
                                        [~,aploc] = max(hev(chtr(1):chtr(2)));
                                        rrp = sort([rrp; tt(round(chtr(1)+aploc-1))-tt(1)]);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                
                
                npk = numel(rrp);
                if npk>0
                    
                    for np = 1:npk
                        
                        if shuffle_peaks
                            % generate random peaks for permutation control
                            % condition
                            
                            sheps = NaN(shuffle_peaks, neegchans, nhepinds);
                            for sp = 1:shuffle_peaks
                                ploc = round(diff(bef_aft)*rand*srate);
                                ff = round(ploc+indrange_hep - pulsedelay);
                                if ff(1)>0 && ff(2)<zind ... % prestimulus
                                        || ff(1)>zind && ff(2)<nerpinds % poststimulus
                                    sheps(sp,:,:) = ep2(selchans, ff(1):ff(2));
                                end
                            end
                            
                            if diff(ff)==nhepinds, ff(2) = ff(2)-1; end
                            if ff(1)>0 && ff(2)<zind
                                bhep = [bhep; nanmean(sheps,1)];
                                hb(1) = hb(1) + 1;
                                
                            elseif ff(1)>zind && ff(2)<nerpinds
                                ahep = [ahep; nanmean(sheps,1)];
                                hb(2) = hb(2) + 1;

                            end
                            
                        else
                            % the actual peaks from the data

                            pkind = round(rrp(np)*srate);
                            if diff(hev([pkind pkind+1]))>0
                                % when the detected peak is slightly before
                                % the detected peak
                                
                                if pkind+pkrange(2)<=nerpinds
                                    % ensure it is less than the span
                                    % within which to calculate the actual
                                    % peak when the peak is slightly
                                    % earlier than the detected peak
                                    
                                    [pk,ploc] = max(hev(pkind:pkind+pkrange(2)));
                                    ploc = ploc + pkind -1;
                                    
                                    if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'ro'), end
                                    
                                    ff = round(ploc+indrange_hep - pulsedelay);
                                    if diff(ff)==nhepinds, ff(2) = ff(2)-1; end
                                    if ff(1)>0 && ff(2)<zind
                                        bhep = [bhep; ep2(selchans, ff(1):ff(2))];
                                        hb(1) = hb(1) + 1;
                                        if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'go'), end
                                        
                                    elseif ff(1)>zind && ff(2)<nerpinds
                                        ahep = [ahep; ep2(selchans, ff(1):ff(2))];
                                        hb(2) = hb(2) + 1;
                                        if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'go'), end
                                    end
                                end
                                
                            else
                                if pkind-pkrange(1)>0
                                    pk = max(hev(pkind-pkrange(1):pkind));
                                    ploc = find(hev==pk);
                                    
                                    if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'ro'), end
                                    
                                    ff = round(ploc+indrange_hep - pulsedelay);
                                    
                                    if ff(1)>0 && ff(2)<zind
                                        bhep = [bhep; ep2(selchans, ff(1):ff(2))];
                                        hb(1) = hb(1) + 1;
                                        if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'go'), end
                                        
                                    elseif ff(1)>zind && ff(2)<nerpinds
                                        ahep = [ahep; ep2(selchans, ff(1):ff(2))];
                                        hb(2) = hb(2) + 1;
                                        if do_plots==1, plot(bef_aft(1)+ploc/srate, pk, 'go'), end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        if ~isempty(bhep)
            bhep = reshape(bhep, [neegchans hb(1) nhepinds]);
            if do_plots==2,plot(squeeze(mean(bhep,1))'), end
            ep1(nt, 1, :, :) = permute(nanmean(bhep, 2), [1 3 2]);
            
            if save_indiv_hep
                % save the individual heps and their corresponding trial
                % numbers
                ep3 = [ep3; permute(bhep, [2 1 3])];
                trialnums = [trialnums repmat(nt, [1 size(bhep,2)])];
            end
        end
        hb_posneg(nt,pn) = hb(1);
        
        if ~isempty(ahep)
            ahep = reshape(ahep, [neegchans hb(2) nhepinds]);
            ep1(nt, 2, :, :) = permute(nanmean(ahep, 2), [1 3 2]);
        end
        
    end
    
    if blsub_flag
        %%% do baseline subtraction
        ep1 = ep1 - repmat(squeeze(mean(ep1(:, :, :, baseinds), 4)), ...
            [1 1 1 nhepinds]);
        
        if save_indiv_hep
            ep3 = ep3 - repmat(squeeze(mean(ep3(:, :, baseinds), 3)), ...
                [1 1 nhepinds]);
        end
    end
    %     val_trials{pn} = vt;
    epochs{pn} = ep1;
    if save_indiv_hep
        indiv_hep{pn} = ep3;
        indiv_trialnums{pn} = trialnums;
    end
end


if ~nargout
    %%% find the channel indices
    nselchan = numel(selchans);
    chansel = false(1, EEG.nbchan);
    chanlabs = {EEG.chanlocs.labels};
    for ns = 1:nselchan
        chansel = chansel|strcmp(chanlabs, selchans{ns});
    end
    
    plotcols = 'br';
    figure, hold on
    for pn = 1:2
        plot(tt, squeeze(mean(mean(epochs{pn}(:, chansel, :), 1), 2)), ...
            plotcols(pn), 'LineWidth', 2)
    end
    legend({'pos', 'neg'})
end
end