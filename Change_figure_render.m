figfiles = dir('*.fig');

for f = 1:1:length(figfiles)
    open(figfiles(f).name)
    F = gcf;
    F.Renderer
    F.Renderer = 'painters';
    savefig(F,figfiles(f).name)
    svgname = strsplit(figfiles(f).name,'.');
    saveas(F,svgname{1},'svg')
    close(F)
end