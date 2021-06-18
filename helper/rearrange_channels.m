function [chansel, chansavail] = rearrange_channels(chanlabs, selchans, ...
    nselchan)
%
% chanlabs: channels of the current dataset
% selchans: set of channels in whose order chanlabs indices will be
% align/rearrange to
% nselchans: first nselchans to keep, remaining to discard

% rearrange channel numbers to ensure order of channels is
% consistent across subjects

chansavail = true(1, nselchan);
chansel = [];
for nsch = 1:nselchan
    findchan = find(strcmp(chanlabs, selchans{nsch}));
    if ~any(findchan)
        chansavail(nsch) = false;
    else
        chansel = [chansel findchan];
    end
end

end