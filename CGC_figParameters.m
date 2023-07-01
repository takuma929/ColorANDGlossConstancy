% script to save parameters for figure generation
function CGC_figParameters(repo_basedir)
    twocolumn = 18.5; % size for for two-column figure (Journal of Vision)
    onecolumn = twocolumn/2; % size for for one-column figure (Journal of Vision)

    fontsize = 7; % fontsize for general use
    fontsize_axis = 8; % fontsize for axis label
    fontname = 'Arial'; % Use Arial for font

    % save figure parameter
    save(fullfile(repo_basedir,'data','CGC_FigParameters'),'twocolumn','onecolumn','fontsize','fontsize_axis',...
        'fontname','repo_basedir')
end