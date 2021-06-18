function hepsret_rt(exp_medriv)

sessdir = exp_medriv.session_dir;
cd(sessdir);

remove_outliers = 1;

for ns = 1:exp_medriv.nsubj
    subj_data = exp_medriv.data(ns);
    bdir = fullfile(sessdir, subj_data.dir_name, 'Behavior', ...
        num2str(subj_data.subj_code));
    if isfolder(bdir)
        cd(bdir);
        sret_img = importdata('SRET_.txt');
        numwords = size(sret_img.data, 1);
        if numwords==60
            neg_words = sret_img.data(:, 2)>30;
            pos_words = sret_img.data(:, 2)<=30;
        else
            disp('err')
        end
        valid_trials = sret_img.data(:, 3)~=9999;
        if sum(neg_words)~=sum(pos_words)
            sprintf('pos and neg words not equal!')
        end
        
        sret_beh=importdata('iPanda_SRET copy.log', '\t', 10); % read the response file
        if size(sret_beh.textdata, 1)<50
            sret_beh=importdata('iPanda_SRET copy.log', '\t', 12); % read the response file
        end
        
        % the actual experiment always starts with the first column being
        % '25', so remove anything before that
        firstline = find(strcmp(sret_beh.textdata(:, 1), '25'));
        sret_beh.textdata = sret_beh.textdata(firstline:end, :);
        
        % find all the rows that have picture-picture in the 2nd and 3rd
        % column - these are when the trial was presented
        p = find((strcmp(sret_beh.textdata(:, 2), 'Picture')) ...
            & (strcmp(sret_beh.textdata(:, 3), 'Picture'))); 

        % find the rows that have picture-picture in the 3rd column in
        % consecutive rows - these are trials when no response was given
        p1 = sret_beh.textdata(:, [2 3 5]);
        pp = strcmp(p1(1:end-1, 2), 'Picture') & strcmp(p1(2:end, 2), 'Picture');
        
        if subj_data.subj_code==4165
            pp(end) = true;
        end
        
        if sum(pp)~=sum(~valid_trials)
            % number of invalid trials should be the same as when there is
            % a picture-picture
            %             fprintf('invalid: trials don''t match up\n')
           
        end
        
        p2 = p1(p+1, [2 3]); % get the one next to picture picture as that is the
        p3 = strcmp(p2(:, 1), '201') | strcmp(p2(:, 1), '202');
        p2 = p2(:, 2);

        resp_time = NaN(1, numwords);
        for np = 1:numel(p3)
            if p3(np)
                resp_time(np) = str2double(p2{np});
                %             else
                %                 fprintf('invalid trial\n')
            end
        end
        
        resp_time = resp_time(valid_trials);        
        pos_words = pos_words(valid_trials);
        neg_words = neg_words(valid_trials);
        
        if remove_outliers
            notout = ~isoutlier(resp_time, 'median');
            resp_time = resp_time(notout);
            pos_words = pos_words(notout);
            neg_words = neg_words(notout);
        end
        posmean(ns) = mean(resp_time(pos_words));
        posmedian(ns) = median(resp_time(pos_words));
        negmean(ns) = mean(resp_time(neg_words));
        negmedian(ns) = median(resp_time(neg_words));
        allposmean(ns) = mean(resp_time);
        
    end
end

median_flag = 0;

if median_flag
    posm = posmedian;
    negm = negmedian;
else
    posm = posmean;
    negm = negmean;
end

allposmean = allposmean(logical(allposmean));
fprintf('\n reaction time across conditions: mean = %g, stderr = %g\n',...
    mean(allposmean), std(allposmean)/sqrt(numel(allposmean)))
end