function [V,Vc,linebyline]=readobj(fname,fpath)
%This simple OBJ parsing routine will work with ASCII encoded files only. 
%It loads vertex coordinates only, and preserves the rest of the file unchanged.
%Vertex colors, which may be appended to the same line as the vertex coordinates,
%are stored in a separate matrix Vc for writing back into the file.
%2D or 4D vertices (x-y-z-w) are not supported. It should work with
%multiple meshes, just as long as all the meshes belong to the same drum
%and are aligned correctly relative to one another.

fid = fopen(fullfile(fpath,fname),'r');
if fid<0
    error(['Cannot open ' fname '.']);
end

l = textscan(fid,'%s','delimiter','\n');
linebyline = l{1};
fprintf('\nRead %d lines from %s', length(linebyline),fname);
fclose(fid);

vertex_lines = linebyline(strncmp('v ',linebyline,2),:);
V = zeros(length(vertex_lines), 3);
[~, lenVs] = sscanf(vertex_lines{1}, '%*c%f'); %Check if there are 3D vertex coordinates, sometimes followed by vertex color RGB values

if lenVs == 3
    Vc = -1;
    for i = 1:length(vertex_lines)
        v = sscanf(char(vertex_lines{i}), 'v %f %f %f');
        V(i, :) = v';
    end
elseif lenVs == 6
    Vc = zeros(length(vertex_lines), 3);
    for i = 1:length(vertex_lines)
        v = sscanf(char(vertex_lines{i}), 'v %f %f %f %f %f %f');
        V(i, :) = v(1:3)';
        Vc(i, :) = v(4:6)';
    end
else
    error('No 3D vertex data encountered.');
end

return