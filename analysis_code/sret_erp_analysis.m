function epochs = sret_erp_analysis(EEG, eeg_mode, posneg_wordinds, bef_aft, ...
    selchans, smooth_erp, blsub_flag)

if ~exist('blsub_flag', 'var')||isempty(blsub_flag), blsub_flag = true; end
if ~exist('smooth_erp', 'var')||isempty(smooth_erp), smooth_erp = 0; end
if ~exist('eeg_mode', 'var')||isempty(eeg_mode), eeg_mode = 1; end
trigsvals = cell(1, 2);
artif_trigs = {'rej1', 'rej2'};

switch eeg_mode
    case 1
        % for actichamp addS to the triggers
        for pn = 1:2
            tinds = posneg_wordinds{pn};
            ntind = numel(tinds);
            winds = cell(1, ntind);
            for nt = 1:ntind
                winds{nt} = add_actichamp_S(tinds(nt));
            end
            trigsvals{pn} = winds;
        end
        
    case 2
        % for emotiv convert directly to strings
        for pn = 1:2
            trigsvals{pn} = num2strcell(posneg_wordinds{pn});
        end
end

%%% get epochs corresponding to pos and neg words
srate = EEG.srate;
nerpinds = round(diff(bef_aft*srate))+1;
baseinds = 1:diff([bef_aft(1) 0]*srate);
epochs = cell(1, 2);
val_trials = cell(1, 2);
for pn = 1:2
    tinds = trigsvals{pn};
    ntind = numel(tinds);
    winds = cell(1, ntind);
    vt = true(1, ntind);
    ep1 = zeros(ntind, EEG.nbchan, nerpinds);
    for nt = 1:ntind
        [ep, ~, artif_mask] = Get_Epochs(EEG, tinds{nt}, [], ...
            floor(bef_aft*srate), [], artif_trigs);
        if ~isempty(ep)
            ep2 = ep{1};
            if smooth_erp
                for nn = 1:size(ep2, 1)
                    ep2(nn, :) = smooth(ep2(nn, :), smooth_erp*srate/1000);
                end
            end
            ep1(nt, :, :) = ep2;
            vt(nt) = sum(artif_mask{1})>size(ep{1}, 2)-10;
        else
            vt(nt) = false;
            sprintf('notrial pn=%g nt=%g', pn, nt)
        end
    end
    
    if blsub_flag
        %%% do baseline subtraction
        ep1 = ep1 - repmat(mean(ep1(:, :, baseinds), 3), [1 1 nerpinds]);
    end
    val_trials{pn} = vt;
    epochs{pn} = ep1;
    trigsvals{pn} = winds;
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
    tt = linspace(bef_aft(1), bef_aft(2), nerpinds);
    for pn = 1:2
        plot(tt, squeeze(mean(mean(epochs{pn}(:, chansel, :), 1), 2)), ...
            plotcols(pn), 'LineWidth', 2)
    end
    legend({'pos', 'neg'})
end
end