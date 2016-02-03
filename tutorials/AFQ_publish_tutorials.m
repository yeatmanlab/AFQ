function AFQ_publish_tutorials(remotename)
% Publish all the AFQ tutorials as html pages push to gh-pages
%
% AFQ_publish_tutorials(remotename)

if ~exist('remotename', 'var') || isempty(remotename)
    remotename = 'origin';
end
afqdir = AFQ_directories;
tdir = fullfile(afqdir,'tutorials');

% List of tutorials to publish
tlist = {'AFQ_example.m')

opts.format = 'html';
opts.outputDir = tdir;
for ii = 1:length(tlist)
    publish(tlist{ii},opts); 
end

% Checkout gh-pages branch
cd(afqdir)
system('git checkout gh-pages')
system('git add tutorials')
system('git commit -m ''changes to tutorials''')
system(sprintf('git push %s gh-pages', remotename))
system('git checkout master')