function exp_hepsret = hepsret_experiment

exp_hepsret.session_dir = '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data';

exp_hepsret.eeg_datadir = 'Data/eeglab';
exp_hepsret.eeg_epochdir = 'Epochs';

subj_codes = [ 4146 4147 4148 4149 4151 4152 ...
    4153 4154 4155 4156 4158 4159 4160 4161 4162  ...
    4164 4165 4167 4168 4169 4170  4173 4181 4186 ...
    4190 4191 4192 4193];

years_practice = [ 46 43 40 17 38 27,...
    0 0 35 15 38 11 0 29 38 , ...
    0 0 0 0 0 28  8 0 0, ...
    0 0 0 0];

nsubj = numel(subj_codes);

subj_age = [ 66 69 64 37 58 54,...
    51 60 65 33 54 43 43 62.5 62  ...
    73 29 57 60 32 62.5  34 62 56 ...
    65 29 52 68];

subj_sex = [ zeros(1, 4) ones(1, 2), ...
    1 0 0 1 0 1 1 0 0 , ...
    0 1 1 0 0 1  0 1 0, ...
    0 1 0 0]; % 0 males, 1 females

subj_edu = [ 20 12 18 20 12 18 ... % check 4145
    18 4154 18 16 18 18 18 14 14  ...
    18 14 18 16 14 16  21 18 20 ...
    14 14 14 16];

subj_valid_sret = 2*ones(1, nsubj);

subj_valid = true(1, nsubj);

exp_hepsret.nsubj = nsubj;

exp_hepsret.hbchanname = 'BVP';
exp_hepsret.hrchanname = 'HRBVP';

cd(exp_hepsret.session_dir)
datadirs = dir;
datadirs = datadirs([datadirs.isdir]);
datadirs = datadirs(3:end);

alldirs = {datadirs.name};

for ii = 1:nsubj
    
    findirs = strfind(alldirs, [num2str(subj_codes(ii)) '_']);
    dirname = alldirs{~isemptycell(findirs)};
    
    exp_hepsret.data(ii).subj_code = subj_codes(ii);
    exp_hepsret.data(ii).subj_age = subj_age(ii);
    exp_hepsret.data(ii).subj_edu = subj_edu(ii);
    exp_hepsret.data(ii).subj_sex = subj_sex(ii);
    exp_hepsret.data(ii).years_practice = years_practice(ii);
    exp_hepsret.data(ii).subj_valid = subj_valid(ii);
    exp_hepsret.data(ii).subj_valid_sret = subj_valid_sret(ii);
    
    exp_hepsret.data(ii).dir_name = dirname;
    sdate = dirname(6:11);
    sdate = sdate([3 4 1 2 5 6]);
    exp_hepsret.data(ii).session_date = sdate;
    
end
end