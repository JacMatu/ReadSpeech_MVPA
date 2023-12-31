% (C) Copyright 2022 Remi Gau

cwd = fileparts(mfilename('fullpath'));
addpath(fullfile(cwd, 'lib', 'bidspm'));
bidspm init;

bids_dir = fullfile(cwd, '..', 'inputs', 'ds003717');
output_dir = fullfile(cwd, '..', 'outputs');

participants = {'01', '02'};

%%
bidspm(bids_dir, output_dir, 'subject', ...
       'action', 'preprocess', ...
       'participant_label', participants, ...
       'dry_run', false, ...
       'verbosity', 2, ...
       'space', {'individual', 'IXI549Space'}, ...
       'task', {'SESS01'}, ...
       'dummy_scans', 0, ...
       'ignore', {}, ...
       'fwhm', 6, ...
       'skip_validation', true);
