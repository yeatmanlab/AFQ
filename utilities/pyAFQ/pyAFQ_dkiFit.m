function pyAFQ_dkiFit(dwi_data, bval_files, bvec_files, mask, out_dir)
% Fit DKI parameters using pyAFQ
%
% pyAFQ_dkiFit(dwi_data, bval_files, bvec_files, mask, out_dir)
%
% This function makes a command line call to pyAFQ_dki to fit the dki model
% and write out dki parameters as nifti images. For this to run you must
% have pyAFQ and dipy installed. After cloning these github repos change 
% into the directorys and:
% python setup.py develop
%Additionally you should have a patht o
% your local bin directory in your .bashrc file
% export PATH="/home/jyeatman/.local/bin/:$PATH"

%% Create strings for the different input files
% We want the files separated by commas but no spaces
d = dwi_data{1}
for ii = 2:length(dwi_data)
    d = horzcat(d, ',',dwi_data{ii}); 
end

l = bval_files{1}
for ii = 2:length(bval_files)
    l = horzcat(l, ',', bval_files{ii}); 
end

c = bvec_files{1}
for ii = 2:length(bvec_files)
    c = horzcat(c, ',' , bvec_files{ii}); 
end

%% Make the system call
cmd = sprintf('pyAFQ_dki -d %s -l %s -c %s -m %s -o %s',d, l, c, mask, out_dir)
system(cmd);
