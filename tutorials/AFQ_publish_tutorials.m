function AFQ_publish_tutorials(remotename)
% Publish all the AFQ tutorials as html pages push to gh-pages
%
% AFQ_publish_tutorials(remotename)
%
% inputs:
% remotename = string for a git remote repository. Default is 'origin'
%
% Written by JDY, Feb, 2016

if ~exist('remotename', 'var') || isempty(remotename)
    remotename = 'origin';
end
afqdir = AFQ_directories;
tdir = fullfile(afqdir,'tutorials');
cd(afqdir)

% List of tutorials to publish
tlist = {'AFQ_example.m', 'AFQ_example_GroupComparison.m'};

% Publish tutorial list
opts.format = 'html';
opts.outputDir = fullfile(tdir,'htmltmp');
for ii = 1:length(tlist)
    publish(tlist{ii},opts); 
end

% Checkout gh-pages branch add tutorials, push to github and change back to
% master branch
[~, curbr] = system('git symbolic-ref --short HEAD')
system('git checkout gh-pages')
system(sprintf('mv %s %s',fullfile(tdir,'htmltmp','*'),tdir))
system(sprintf('rm -d %s %s',fullfile(tdir,'htmltmp')))
system('git add tutorials')
system('git commit -m ''changes to tutorials''')
system(sprintf('git push %s gh-pages', remotename))
system(sprintf('git checkout %s',curbr))