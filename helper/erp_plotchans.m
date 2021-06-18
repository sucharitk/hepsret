function erp_plotchans(erp_data, chanlabs, timeval, ax,...
    within_err, separate_c1)
%
% function erp_plotchans(erp_data, chanlabs, timeval)
%
% plots erps for each channel on different subplots
%
% erp_data = subj x channels x time points [x condition 1 x condition 2]
% condition 1 is subplotted on the same figure
% condition 2 is plotted on a new figure
% chanlabs: channel labels
% timeval: time points
% ax     : specify axis
% within_err: plot within c1 condition errorbars (loftus and mason)
% separate_c1 : plot yes and no on separate figures for ungrouping


nd = ndims(erp_data);
szd = size(erp_data);
if nd>4, ncd2 = szd(5); else, ncd2=1; end
if nd>3, ncd1 = szd(4); else, ncd1=1; end
nsubj = szd(1); nchans = szd(2); ntp = szd(3);
if ~exist('timeval', 'var'), timeval = 1:ntp; end % hack the time to inds if timepoint values are not specified
if ~exist('within_err', 'var'), within_err = false; end % hack the time to inds if timepoint values are not specified
if ~exist('separate_c1', 'var'), separate_c1 = false; end % hack the time to inds if timepoint values are not specified


nsubs = ceil(sqrt(nchans)); % num of subplots
if nsubs*(nsubs-1) >= nchans
    nsubs = [nsubs-1 nsubs];
else
    nsubs = [nsubs nsubs];
end

plotcols = {'b', 'r', 'g'; 'b-.', 'r-.', 'g'};


for n2 = 1:ncd2
    if ~separate_c1
        figure('Name', sprintf('c2 = %g', n2))
    end
    for nc = 1:nchans
        subplot(nsubs(1), nsubs(2), nc)
        hold on
        if within_err
            % calculate within subject error between conditions
            e2 = squeeze(diff(erp_data(:,nc,:,[2 1],n2), [], 4));
            smm = nanstd(e2,1)/sqrt(nsubj);
        end
        for n1 = 1:ncd1
            if separate_c1
                figure('Name', sprintf('c1 = %g, c2 = %g', n1, n2))
            end

            mm = squeeze(erp_data(:, nc, :, n1, n2));
            mmm = nanmean(mm, 1);
            %             plot(timeval, mmm, plotcols{1, n1}, 'LineWidth', 2)
            if ~within_err
                smm = nanstd(mm, 1)/sqrt(nsubj);
            end
            shadedErrorBar(timeval, mmm, smm, plotcols{1, n1},0)

            if exist('ax','var'),axis(ax);end
        end
        title(chanlabs{nc})
    end
end

end