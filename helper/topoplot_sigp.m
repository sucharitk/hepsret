function topoplot_sigp(data,chanlocs,sigchans,maplims,labelstring,...
    scale_mark,labelposoffset,cbar, grayscale)

if ~exist('labelstring','var'), labelstring = ''; end
if ~exist('scale_mark','var'), scale_mark = 1; end
if ~exist('cbar','var'), cbar = false; end
if ~exist('grayscale','var'), grayscale = false; end
if ~exist('maplims','var'), maplims = [min(data(:)) max(data(:))]; end

topoplot(data, chanlocs, ...
    'maplimits', maplims , 'emarker2', {sigchans,'*','w', 18*scale_mark, 2});
topoplot(data, chanlocs, ...
    'maplimits', maplims , 'emarker2', {sigchans,'o','k', 17*scale_mark, 2});
topoplot(data, chanlocs, ...
    'maplimits', maplims , 'emarker2', {sigchans,'.','k', 10*scale_mark});
if cbar
    c = colorbar('Ticks', [-0 4], 'FontSize', 20);
    c.Label.String = labelstring;
    % c.Label.FontSize = 30;
    % c.Label.Rotation = 270;
    c.Label.Position(1) = c.Label.Position(1)-labelposoffset;
    % x1=get(gca,'position');
    % x=get(c,'position');
    % x([2 3 4]) = [.19 0.035 .65];
    % set(c, 'Position', x)
    % set(gca, 'position', x1)
end
if grayscale
    colormap gray
end

end