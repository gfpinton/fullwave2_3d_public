function [coords] = mapToCoords3D (map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebecca Jones
% convert 3D map to coordinate vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X Y Z]=ind2sub(size(map),find(map~=0));
coords = [X Y Z];
