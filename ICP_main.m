%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main script for executing the inside-ICP fitting of column
% drums in Matlab. Instructions are within the code.
% Authors Philip Sapirstein and Eric Psota
% University of Nebraska-Lincoln, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;

%First, place the scan files to be analyzed into a subfolder from where
%this script is located. Files will be written out to the current path.
scan_dir = 'DrumScans';

%Three models are supplied with the code. Only one can be selected, so only
%drums with similar flute depths can analyzed at once. These files must be
%located in the same directory as the scripts.
model_name = '_Idealized_drum-RegularFlutes.obj';
%model_name = '_Idealized_drum-DeepFlutes.obj';
%model_name = '_Idealized_drum-Faceted.obj';

[script_path,~,~] = fileparts(mfilename('fullpath'));
objfilelist = dir(fullfile(script_path, scan_dir, '*.obj'));
scan_names = {objfilelist.name}';
if isempty(scan_names)
    error('Invalid directory, or no .OBJ files present.');
end

%Name the file for saving the estimated parameters here; it will be
%appended if it already exists.
fid = fopen('_ModelFittingICP.csv','at+');
fprintf(fid, 'Process initiated %s\n', datestr(now));
fprintf(fid, 'Model Name, Height, Top Radius, Mid Radius, Bottom Radius, Taper Angle, Err, MaxH\n');

paralleljob = gcp('nocreate');
if isempty(paralleljob)
    paralleljob = parpool('local8',6);
elseif paralleljob.NumWorkers < 6
    delete(paralleljob);
    paralleljob = parpool('local8',6);
end

for idx = 1:size(scan_names,1)
    %All the work is done by the following subroutine, which also saves
    %figures and a reoriented copy of the scan in the active folder
    [Parameters] = ICPFitDrumToModel(model_name, scan_names{idx}, fullfile(script_path, scan_dir));
    
    %Save the estimated parameters for each scan
    fprintf(fid, '%s, %f, %f, %f, %f, %f, %f, %f\n', scan_names{idx}, Parameters(1), Parameters(2), Parameters(3), Parameters(4), Parameters(5), Parameters(6), Parameters(7));
end
fprintf(fid, 'Process concluded %s\n', datestr(now));
fclose(fid);
delete(paralleljob);
