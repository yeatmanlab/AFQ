function test_pyAFQ
% Test AFQ calls to pyAFQ for compatibility
%
% test_pyAFQ
%
% Concurent development is being done in AFQ (matlab) and pyAFQ (python).
% These two packages should be fully compatible. This function will test
% that AFQ is compatible with the current version of pyAFQ.

%% Test dki fitting implimented in pyAFQ

% Process test data with pyAFQ_dkiFit
d = fullfile(which('test_pyAFQ', 'data'));
dwi_data   = fullfile(d, 'dwi_testdki.nii.gz');
bval_files = fullfile(d, 'dwi_testdki.nii.gz');
bvec_files = fullfile(d, 'dwi_testdki.nii.gz');
mask       = fullfile(d, 'dwi_mask.nii.gz');
out_dir     = d;
try
    pyAFQ_dkiFit(dwi_data, bval_files, bvec_files, mask, out_dir);
catch
    error('\n Call to pyAFQ_dki failed\n');
end

% Compare output files to expectation

