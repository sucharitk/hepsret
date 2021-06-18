function hepsret_figures(exp_hepsret, stats, fig_num, topo_time)
%
% 

cd(exp_hepsret.session_dir)
load('chanlocs.mat')
switch fig_num
        
    case 1
        % plot the 3 panels of Fsigure 3 for the paper
        load('hepsret_datafiles/hrep.mat')
        
        data(1, :, :, :) = squeeze(hrep(:, 2, 1, :, :));
        data(2, :, :, :) = squeeze(hrep(:, 2, 2, :, :));
        nht = size(hrep,5);
        tt = linspace(trange_hep(1), trange_hep(2), nht);
        
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        
        if exist('topo_time', 'var')&&~isempty(topo_time)
            % if there is a specified time over which to calculate the
            % topographies and the ERPs
            timestat = stats.time;
            
            sigtimes = topo_time;
            sigchans = all(sigclust(:,timestat>=sigtimes(1) &...
                timestat<=sigtimes(2)),2);
        else
            % otherwise just do an OR and average over the channels and
            % times that are significant
            sigtimes = stats.time(logical(sum(sigclust,1)));
            sigchans = logical(sum(sigclust,2));
        end
        
        topodata = squeeze(mean(diff(mean(data([2 1],:,:,...
            tt>=sigtimes(1) & tt<=sigtimes(2)),2),[],1), 4));
        chans2plot = get_channels_from_labels({chanlocs.labels}, labchans);
        figure, topoplot_sigp(topodata', chanlocs(chans2plot),...
            find(sigchans))

        within_err = true;
        seprate_c1 = true; % plot yes and no on separate figures for ungrouping
        erp_plotchans(permute(mean(hrep(:,:,:,sigchans,:),4), [1 4 5 3 2]), ...
            labchans(1), tt, [-.1 .5 -2 7],...
            within_err, seprate_c1);
        

    case 2
        % control analysis - shuffled peaks
        load('hepsret_datafiles/hrep_shuffle.mat')
        
        nht = size(hrep,5);
        tt = linspace(trange_hep(1), trange_hep(2), nht);
        
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        sigchans = logical(sum(sigclust,2));
            
        erp_plotchans(permute(mean(hrep(:,:,:,sigchans,:),4), [1 4 5 3 2]), ...
            labchans(1), tt, [-.1 .5 -2 7]);
end