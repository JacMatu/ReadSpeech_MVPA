function value = moae_dir()
  value = bids.internal.file_utils(fullfile(get_test_data_dir(), ...
                                            '..', '..', ...
                                            'demos', 'spm', 'moae'), ...
                                   'cpath');
end
