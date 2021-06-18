function [tstat,pval,diffval]=hepsret_t_hepdur(exp_hepsret, stats, ...
    stat_type, reg_or_shuf)
%
% do ttest in the channelxtime cluster for different durations of hep
% calculation interval
%

if ~exist('reg_or_shuf', 'var'), reg_or_shuf = 1; end

cd(exp_hepsret.session_dir)
switch reg_or_shuf
    case 1
        load('hepsret_datafiles/hrep.mat')
    case 2
        load('hepsret_datafiles/hrep_shuffle.mat')
end

switch stat_type
    case 1
        data(1, :, :, :) = squeeze(hrep(:, 1, 1, :, :));
        data(2, :, :, :) = squeeze(hrep(:, 1, 2, :, :));
        nht = size(hrep,5);
        tt = linspace(trange_hep(1), trange_hep(2), nht);
        
        data = data(:,:,:,tt>=0);
        
        
        nsubj = size(data,2);
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        
        d2 = NaN(2, nsubj);
        for yn = 1:2
            for ns = 1:nsubj
                d1 = squeeze(data(yn,ns,:,:));
                d2(yn,ns) = mean(d1(sigclust));
            end
        end
        
        [~,pval,~,stat]=ttest(diff(d2));
        tstat = stat.tstat;
        fprintf('t = %g, p = %g\n', tstat, pval)
        
    case 2
        % unpleasant HEPs
        data(1, :, :, :) = squeeze(hrep(:, 2, 1, :, :));
        data(2, :, :, :) = squeeze(hrep(:, 2, 2, :, :));
        nht = size(hrep,5);
        tt = linspace(trange_hep(1), trange_hep(2), nht);
        
        data = data(:,:,:,tt>=0);
        
        nsubj = size(data,2);
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        
        d2 = NaN(2, nsubj);
        for yn = 1:2
            for ns = 1:nsubj
                d1 = squeeze(data(yn,ns,:,:));
                d2(yn,ns) = mean(d1(sigclust));
            end
        end
        
        diffval = diff(d2([2 1], :));
        [~,pval,~,stat]=ttest(diffval);
        tstat = stat.tstat;
        %         diffval = mean(dd2);
        fprintf('t = %g, p = %g\n', tstat, pval)
        
    case 3
        % pos heart rate
        
        data(1, :, :) = squeeze(hr_sret(:, 1, 1, :));
        data(2, :, :) = squeeze(hr_sret(:, 1, 2, :));
        nht = size(hr_sret,4);
        tt = linspace(bef_aft_hr(1), bef_aft_hr(2), nht);
        
        data = data(:,:,tt<0);
        
        
        nsubj = size(data,2);
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        
        d2 = NaN(2, nsubj);
        for yn = 1:2
            for ns = 1:nsubj
                d1 = squeeze(data(yn,ns,:,:));
                d2(yn,ns) = mean(d1(sigclust));
            end
        end
        
        [~,pval,~,stat]=ttest(diff(d2));
        tstat = stat.tstat;
        fprintf('t = %g, p = %g\n', tstat, pval)
        
    case 4
        % negheart rate
        
        data(1, :, :) = squeeze(hr_sret(:, 2, 1, :));
        data(2, :, :) = squeeze(hr_sret(:, 2, 2, :));
        nht = size(hr_sret,4);
        tt = linspace(bef_aft_hr(1), bef_aft_hr(2), nht);
        data = data(:,:,tt<0);
        
        [~,pval,~, stat]=ttest(diff(mean(data,3)));
        tstat = stat.tstat;
        fprintf('HR diff neg words: t = %g, p = %g\n', tstat, pval)
        
    case 2
        % unpleasant HEPs
        data(1, :, :, :) = squeeze(hrep(:, 2, 1, :, :));
        data(2, :, :, :) = squeeze(hrep(:, 2, 2, :, :));
        nht = size(hrep,5);
        tt = linspace(trange_hep(1), trange_hep(2), nht);
        
        data = data(:,:,:,tt>=0);
        
        nsubj = size(data,2);
        clustnum = 1;
        sigclust = stats.posclusterslabelmat==clustnum;
        
        d2 = NaN(2, nsubj);
        for yn = 1:2
            for ns = 1:nsubj
                d1 = squeeze(data(yn,ns,:,:));
                d2(yn,ns) = mean(d1(sigclust));
            end
        end
        
        diffval = diff(d2);
        [~,pval,~,stat]=ttest(diffval);
        tstat = stat.tstat;
        %         diffval = mean(dd2);
        fprintf('t = %g, p = %g\n', tstat, pval)
        
                
end
end