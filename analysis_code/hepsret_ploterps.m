function hepsret_ploterps(exp_hepsret, plot_type)

cd(fullfile(exp_hepsret.session_dir))

load('chanlocs32.mat')
switch plot_type
    case {1,2,3,5}
        load('hepsret_datafiles/hrep.mat')
    case 4
        load('hepsret_datafiles/hrep_shuffle.mat')
end

switch plot_type
    case 1
        thep = linspace(trange_hep(1), trange_hep(2), size(hrep,5));
        erp_plotchans(permute(hrep, [1 4 5 3 2]), labchans, thep);
        
    case 2
        erp_data = permute(erp_sret, [1 4 5 3 2]);
        terp = linspace(bef_aft_erp(1), bef_aft_erp(2), size(erp_sret,5));
        erp_plotchans(erp_data, labchans, terp);
        
    case 3
        terp = linspace(bef_aft_hr(1), bef_aft_hr(2), size(hr_sret,4));
        hr_sret = permute(hr_sret, [1 5 4 3 2]); % add an extra channel dimension to use with the erp plot function
        erp_plotchans(hr_sret, {'HR'}, terp);
        
    case 4
        % control analysis - shuffled peaks
        thep = linspace(trange_hep(1), trange_hep(2), size(hrep,5));
        erp_plotchans(permute(hrep, [1 4 5 3 2]), labchans, thep);

    case 5
        % post heps
        thep = linspace(trange_hep(1), trange_hep(2), size(hrep,5));
        erp_plotchans(permute(hrep_post, [1 4 5 3 2]), labchans, thep);

end

end