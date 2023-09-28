% (C) Copyright 2021 CPP ROI developers

function test_suite = test_getAtlasAndLut() %#ok<*STOUT>
  try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions = localfunctions(); %#ok<*NASGU>
  catch % no problem; early Matlab versions can use initTestSuite fine
  end
  initTestSuite;
end

function test_lut_wang()

  [atlasFile, lut] = getAtlasAndLut('wang');

  rmRetinoAtlas();

end

function test_lut_neuromorphometrics()

  [atlasFile, lut] = getAtlasAndLut('neuromorphometrics');

end
