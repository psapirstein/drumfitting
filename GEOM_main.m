%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main script for executing the geometric fitting of column
% drums in Matlab. Instructions are within the code.
% Author Philip Sapirstein, with contributions by Eric Psota
% University of Nebraska-Lincoln, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;

%First, place the scan files to be analyzed into a subfolder from where
%this script is located. Files will be written out to the current path.
scan_dir = 'DrumScans';

%Three models are supplied with the code. They are not used during the
%fitting process, but rather to create a figure comparing the reoriented
%scan to the ideal model scaled to the estimated size of the piece. There
%are three models here to accommodate the different flute depths in the
%material analyzed from Olympia. These files must be located in the same
%directory as the scripts.
model_reg = '_Idealized_drum-RegularFlutes.obj';
model_deep = '_Idealized_drum-DeepFlutes.obj';
model_flat = '_Idealized_drum-Faceted.obj';

[script_path,~,~] = fileparts(mfilename('fullpath'));
objfilelist = dir(fullfile(script_path, scan_dir, '*.obj'));
scan_names = {objfilelist.name}';
if isempty(scan_names)
    error('Invalid directory, or no .OBJ files present.');
end

%Name the file for saving the estimated parameters here; it will be
%appended if it already exists.
fid = fopen('_ModelFittingGeom.csv','at+');
fprintf(fid, 'Process initiated %s\n', datestr(now));
fprintf(fid, 'Scan Name, Height, Top Radius, Mid Radius, Bottom Radius, Taper Angle, Flute Depth, E\n');

nFlutes = 20; %Assume 20 flutes, though this can be altered for nonstandard flute counts (e.g. 16)

for idx = 1:size(scan_names,1)
    %All the work is done by the following subroutine, which also saves
    %figures and a reoriented copy of the scan in the active folder
    [Parameters] = GEOMFitDrumByRadius(model_reg, model_deep, model_flat, scan_names{idx}, fullfile(script_path, scan_dir), nFlutes);
    
    %Save the estimated parameters for each scan
    fprintf(fid, '%s, %f, %f, %f, %f, %f, %f, %f\n', scan_names{idx}, Parameters(1), Parameters(2), Parameters(3), Parameters(4), Parameters(5), Parameters(6), Parameters(7));
end
fprintf(fid, 'Process concluded %s\n', datestr(now));
fclose(fid);
