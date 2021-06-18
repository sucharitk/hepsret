function [epochs, latencies, artif_mask] = Get_Epochs(EEG, trig1, trig2, epoch_duration_inds, ...
    upto_next_event, artif_trigs)
%
% [epochs, latencies, artif_mask] = Get_Epochs(EEG, trig1, trig2, epoch_duration_inds, ...
%     upto_next_event, artif_trigs)
%
%
% Reads eeglab structure and gets events ranging between two trigger codes
%
% epoch_duration_inds: (duration of epoch)*srate, if provided then it will
% ignore trig2

if ~exist('upto_next_event', 'var') || isempty(upto_next_event)
    upto_next_event = false; end

evts = {EEG.event.type};
        
if ~isempty(trig1)
    if ischar(trig1)
        calctype = 1;
    else
        if iscellnumeric(evts)
            evts = cell2mat(evts);
            calctype = 2;
        elseif iscellstr(evts)
            calctype = 3;
        else
            calctype = 2;
        end
    end
    
    switch calctype
        case 1
            trig1_evts = find(strcmp(evts, trig1));
        case 2
            trig1_evts = find(evts==trig1);
        case 3
            trig1_evts = find(strcmp(evts, num2str(trig1)));
    end
    numev = numel(trig1_evts);
else
    % If trig1 is empty use trig2 and calculate epochs from the ending to
    % beginning
    if ischar(trig2)
        trig2_evts = find(strcmp(evts, trig2));
    else
        if iscellnumeric(evts), evts = cell2mat(evts); end
        trig2_evts = find(evts==trig2);
    end
    numev = numel(trig2_evts);
end

if ~numev
    % no such event found
    epochs = {};
    latencies = 0;
    artif_mask = 0;
end
rem_artif_epochs = false;

if exist('artif_trigs', 'var')
    if ischar(artif_trigs{1})
        artif_evts{1} = find(strcmp(evts, artif_trigs{1}));
    else
        artif_evts{1} = find(evts==artif_trigs{1});
    end
    if ischar(artif_trigs{2})
        artif_evts{2} = find(strcmp(evts, artif_trigs{2}));
    else
        artif_evts{2} = find(evts==artif_trigs{2});
    end
    if numel(artif_evts{1})~=numel(artif_evts{1})
        % error
        disp('Number of artifact events for the two triggers are not equal')
        epochs = [];
        return
    end
    
    artif_lats{1} = [EEG.event(artif_evts{1}).latency];
    artif_lats{2} = [EEG.event(artif_evts{2}).latency];
    
end

if exist('epoch_duration_inds', 'var') && ~isempty(epoch_duration_inds)
    if numel(epoch_duration_inds)==1, ...
            epoch_duration_inds = [0 epoch_duration_inds-1]; end
    if exist('trig1_evts', 'var')
        % Chop events fromt the start
        epochs = cell(1, numev);
        latencies = zeros(2, numev);
        
        for ii = 1:numev
            t1 = trig1_evts(ii);
            latencies(:, ii) = [round(EEG.event(t1).latency)+epoch_duration_inds(1), ...
                round(EEG.event(t1).latency)+epoch_duration_inds(2)];
            
            
            if exist('artif_trigs', 'var')
                ep = [];
                ae1 = artif_lats{1}(artif_lats{1}>latencies(1, ii) & artif_lats{1}<latencies(2, ii));
                ae2 = artif_lats{2}(artif_lats{2}>latencies(1, ii) & artif_lats{2}<latencies(2, ii));
                nae1 = numel(ae1); nae2 = numel(ae2);

                if  nae1==nae2
                    if ~nae1 || (ae1(1) < ae2(1))
                        ae2 = [latencies(1, ii) ae2];
                        ae1 = [ae1 latencies(2, ii)];
                    end
                elseif nae1 == nae2+1
                    ae2 = [latencies(1, ii) ae2];
                elseif nae2 == nae1+1
                    ae1 = [ae1 latencies(2, ii)];
                else
                    % error
                    error('Number of artifact events for the two triggers are not equal between %d to %d'...
                        , t1, t2)
                    
                end
                
                nae1 = numel(ae1);
                %                 latencies(1, ii) = round(ae2(1));
                %                 latencies(2, ii) = round(ae1(end));
                if rem_artif_epochs
                    for na = 1:nae1
                        ep = [ep EEG.data(:, round(EEG.event(ae2(na)).latency):...
                            round(EEG.event(ae1(na)).latency))];
                        
                    end
                    epochs{ii} = ep;
                else
                    epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
                    artifmask = false(1, size(epochs{ii}, 2));
                    for na = 1:nae1
                        artifmask(round(ae2(na))-latencies(1, ii)+1:...
                            round(ae1(na))-latencies(1, ii)+1) = true;
                    end
                    artif_mask{ii} = artifmask;
                end
                
            else
                epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
            end
            
        end
    else
        % Events from the end
        % Chop events fromt the start
        epochs = cell(1, numev);
        latencies = zeros(2, numev);
        for ii = 1:numev
            t2 = trig2_evts(ii);
            latencies(:, ii) = [round(EEG.event(t2).latency)-epoch_duration_inds(1), ...
                round(EEG.event(t2).latency)+epoch_duration_inds(2)];
            epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
        end
    end
    
else
    
    %     trig2_evts = find(strcmp(evts, trig2));
    switch calctype
        case 1
            trig2_evts = find(strcmp(evts, trig2));
        case 2
            trig2_evts = find(evts==trig2);
        case 3
            trig2_evts = find(strcmp(evts, num2str(trig2)));
    end
    numev2 = numel(trig2_evts);
    
    
    if numev~=numel(trig2_evts)
        if upto_next_event
            % Get events for which we have trig1 upto the next trig2 event
            epochs = cell(1, numev);
            latencies = zeros(2, numev);
                        
            for ii = 1:numev
                t1 = trig1_evts(ii);
                
                if ~isempty(trig2)
                    t2 = trig2_evts(find((trig2_evts-t1)>0, 1)); % Get the very first next event that matches
                else
                    % evaluate to any next trigger doesn't matter if there
                    % is t2
                    t2 = t1+1;
                end
                    
                latencies(:, ii) = [round(EEG.event(t1).latency), ...
                    round(EEG.event(t2).latency)];
                epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
            end
        else
            disp('Number of events for the two triggers are not equal')
            return
        end
    else
        % If equal number of events
        epochs = cell(1, numev);
        latencies = zeros(2, numev);
        for ii = 1:numev
            t1 = trig1_evts(ii);
            t2 = trig2_evts(ii);
            
            if exist('artif_trigs', 'var')
                ep = [];
                ae1 = artif_evts{1}(artif_evts{1}>t1 & artif_evts{1}<t2);
                ae2 = artif_evts{2}(artif_evts{2}>t1 & artif_evts{2}<t2);
                nae1 = numel(ae1); nae2 = numel(ae2);
                if  nae1==nae2
                    if ~nae1 || (ae1(1) < ae2(1))
                        ae2 = [t1 ae2];
                        ae1 = [ae1 t2];
                    end
                elseif nae1 == nae2+1
                    ae2 = [t1 ae2];
                elseif nae2 == nae1+1
                    ae1 = [ae1 t2];
                else
                    % error
                    error('Number of artifact events for the two triggers are not equal between %d to %d'...
                        , t1, t2)

                end
                nae1 = numel(ae1);
                %                 latencies(1, ii) = round(EEG.event(ae2(1)).latency);
                %                 latencies(2, ii) = round(EEG.event(ae1(end)).latency);
                latencies(1, ii) = round(EEG.event(t1).latency);
                latencies(2, ii) = round(EEG.event(t2).latency);
                if rem_artif_epochs
                    for na = 1:nae1
                        ep = [ep EEG.data(:, round(EEG.event(ae2(na)).latency):...
                            round(EEG.event(ae1(na)).latency))];
                        
                    end
                    epochs{ii} = ep;
                else
                    epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
                    artifmask = false(1, size(epochs{ii}, 2));
                    for na = 1:nae1
                        artifmask(round(EEG.event(ae2(na)).latency)-latencies(1, ii)+1:...
                            round(EEG.event(ae1(na)).latency)-latencies(1, ii)+1) = true;
                    end
                    artif_mask{ii} = artifmask;
                end
                
            else
                latencies(:, ii) = [round(EEG.event(t1).latency), ...
                    round(EEG.event(t2).latency)];
                epochs{ii} = EEG.data(:, latencies(1, ii):latencies(2, ii));
            end
        end
    end
    
end
