function [map2] = extendMap3d(map,nbdy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2021-02-28
% LAST MODIFIED: 2021-02-28
% Extend map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map2 = zeros(size(map,1)+2*nbdy,size(map,2)+2*nbdy,size(map,3)+2*nbdy);
% center 
map2(nbdy+1:end-nbdy,nbdy+1:end-nbdy,nbdy+1:end-nbdy) = map;

% faces
for i=1:size(map,1)
  for j=1:size(map,2)
    for k=1:nbdy
      map2(i+nbdy,j+nbdy,k)=map(i,j,1);
    end
    for k=size(map,3)+1+nbdy:size(map,3)+2*nbdy
      map2(i+nbdy,j+nbdy,k)=map(i,j,end);
    end
  end
end

for i=1:size(map,1)
  for k=1:size(map,3)
    for j=1:nbdy
      map2(i+nbdy,j,k+nbdy)=map(i,1,k);
    end
    for j=size(map,2)+1+nbdy:size(map,2)+2*nbdy
      map2(i+nbdy,j,k+nbdy)=map(i,end,k);
    end
    
  end
end

for j=1:size(map,2)
  for k=1:size(map,3)
    for i=1:nbdy
      map2(i,j+nbdy,k+nbdy)=map(1,j,k);
    end
    for i=size(map,1)+1+nbdy:size(map,1)+2*nbdy
      map2(i,j+nbdy,k+nbdy)=map(end,j,k);
    end
    
  end
end


% edges
for i=1:size(map,1)
  map2(i+nbdy,1:nbdy,1:nbdy) = ones(nbdy,nbdy)*map(i,1,1);
  map2(i+nbdy,1:nbdy,end-nbdy+1:end) = ones(nbdy,nbdy)*map(i,1,end);
  map2(i+nbdy,end-nbdy+1:end,1:nbdy) = ones(nbdy,nbdy)*map(i,end,1);
  map2(i+nbdy,end-nbdy+1:end,end-nbdy+1:end) = ones(nbdy,nbdy)*map(i,end,end);
end
for j=1:size(map,2)
  map2(1:nbdy,j+nbdy,1:nbdy)=ones(nbdy,nbdy)*map(1,j,1);
  map2(1:nbdy,j+nbdy,end-nbdy+1:end)=ones(nbdy,nbdy)*map(1,j,end);
  map2(end-nbdy+1:end,j+nbdy,1:nbdy)=ones(nbdy,nbdy)*map(end,j,1);
  map2(end-nbdy+1:end,j+nbdy,end-nbdy+1:end)=ones(nbdy,nbdy)*map(end,j,end);
end
for k=1:size(map,3)
  map2(1:nbdy,1:nbdy,k+nbdy)=ones(nbdy,nbdy)*map(1,1,k);
  map2(1:nbdy,end-nbdy+1:end,k+nbdy)=ones(nbdy,nbdy)*map(1,end,k);
  map2(end-nbdy+1:end,1:nbdy,k+nbdy)=ones(nbdy,nbdy)*map(end,1,k);
  map2(end-nbdy+1:end,end-nbdy+1:end,k+nbdy)=ones(nbdy,nbdy)*map(end,end,k);
end
				% corners
map2(1:nbdy,1:nbdy,1:nbdy)=map(1,1,1);
map2(1:nbdy,1:nbdy,end-nbdy+1:end)=map(1,1,end);
map2(1:nbdy,end-nbdy+1:end,1:nbdy)=map(1,end,1);
map2(1:nbdy,end-nbdy+1:end,end-nbdy+1:end)=map(1,end,end);
map2(end-nbdy+1:end,1:nbdy,1:nbdy)=map(end,1,1);
map2(end-nbdy+1:end,1:nbdy,end-nbdy+1:end)=map(end,1,end);
map2(end-nbdy+1:end,end-nbdy+1:end,1:nbdy)=map(end,end,1);
map2(end-nbdy+1:end,end-nbdy+1:end,end-nbdy+1:end)=map(end,end,end);

