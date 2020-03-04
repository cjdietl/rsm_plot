% Here, the spec file is 'CRO_CD409_RSM.spc' and the myt-files reside in '.\CRO_CD409_RSM\'
% The RSM is associated with scan numer #27. The detector sample distance
% is 228 mm. Discard myt files with a maximum count below 10 to reduce
% data.
out=rsm_map_fourc('CRO_CD409_RSM.spc','.\CRO_CD409_RSM\',27,228,10);

% Use a mesh-interpolation of 1 (amount of interpolated 'cornering' points around an
% actual data point). Project the points along the (010) and (100)
% direction. Use 5E-4 and 1E-3 rsu as a stepsize for the interpolation and
% select a box with the coords (1.9,3.93,2.1,4.1) from this map for
% plotting.

create_gridded_rsm_view(out,1,[0,1,0],[1,0,0],5E-4,1E-3,1.9,3.93,2.1,4.1)