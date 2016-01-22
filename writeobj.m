function writeobj(modelname,V,Vc,linebyline)
%This will save a copy to an ASCII formatted file with new vertex coordinates
%(assumes a 6-decimal precision for the vertex and color values).
fid = fopen(modelname,'w');
vidx = 1;
Vs = cell(size(V,1),1);

if Vc == -1
    for i = 1:size(V,1)
        Vs{i} = sprintf('v %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3));
    end
else
    for i = 1:size(V,1)
        Vs{i} = sprintf('v %.6f %.6f %.6f %.6f %.6f %.6f', V(i,1), V(i,2), V(i,3), Vc(i,1), Vc(i,2), Vc(i,3));
    end
end

for i = 1:size(linebyline,1)
    if strncmp('v ',linebyline{i},2)
        fprintf(fid,'%s\n',Vs{vidx});
        vidx = vidx+1;
    else
        fprintf(fid,'%s\n',linebyline{i});
    end
end

fclose(fid);
end