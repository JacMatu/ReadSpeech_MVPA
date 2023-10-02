function isAKnownAtlas(atlasName)
  % (C) Copyright 2021 CPP ROI developers
  if ~ismember(lower(atlasName), {'anatomy_toobox', ...
                                  'glasser', ...
                                  'hcpex', ...
                                  'neuromorphometrics', ...
                                  'visfatlas', ...
                                  'wang'})
    % TODO throw a proper error here
    error('unknown atlas type');
  end
end