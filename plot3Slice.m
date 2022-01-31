function [] = plot3Slice(mat,idx,varargin)
%plot3Slice(mat,idx,idvec)
%plot3Slice(mat,idx,idvec,coords)

optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin == 1
  idvec = [round(size(mat,1)/2) round(size(mat,2)/2) round(size(mat,3)/2)];
end
if stdargin == 2
  figure(idx)
  idvec = [round(size(mat,1)/2) round(size(mat,2)/2) round(size(mat,3)/2)];
end
if optargin == 1
  idvec = varargin{1};
end

if(length(size(mat)==3))
subplot(2,2,1)
imagesc(squeeze(mat(:,:,idvec(3)))), colorbar, xlabel('2'),ylabel('1')
subplot(2,2,2)
imagesc(squeeze(mat(:,idvec(2),:))), colorbar, xlabel('3'),ylabel('1')
subplot(2,2,3)
imagesc(squeeze(mat(idvec(1),:,:))), colorbar, xlabel('3'),ylabel('2')
end
if(length(size(mat==4)))
  subplot(2,2,1)
  imagesc(squeeze(mat(:,:,idvec(3),:))), colorbar, xlabel('2'),ylabel('1')
  subplot(2,2,2)
  imagesc(squeeze(mat(:,idvec(2),:,:))), colorbar, xlabel('3'),ylabel('1')
  subplot(2,2,3)
  imagesc(squeeze(mat(idvec(1),:,:,:))), colorbar, xlabel('3'),ylabel('2')
end

if optargin == 2
  idvec = varargin{1};
  coords = varargin{2};

  subplot(2,2,1)
  cslice = squeeze(mat(:,:,idvec(3)));
  maxval = max(max(cslice));
  idc = find(coords(:,3)==idvec(3));
  for i=1:length(idc)
    cslice(coords(idc(i),1),coords(idc(i),2)) = maxval;
  end
  imagesc(cslice), colorbar
  axis equal
  xlabel('2'),ylabel('1')
  
  subplot(2,2,2)
  cslice = squeeze(mat(:,idvec(2),:)); % nY and nZ
  idc = find(coords(:,2)==idvec(2));
  maxval = max(max(cslice));
  for i=1:length(idc)
    cslice(coords(idc(i),1),coords(idc(i),3)) = maxval;
  end
  imagesc(cslice), colorbar
  axis equal
  xlabel('3'),ylabel('1')

  subplot(2,2,3)
  cslice = squeeze(mat(idvec(1),:,:)); % nY and nZ
  idc = find(coords(:,1)==idvec(1));
  maxval = max(max(cslice));
  for i=1:length(idc)
    cslice(coords(idc(i),2),coords(idc(i),3)) = maxval;
  end
  imagesc(cslice), colorbar
  axis equal
  xlabel('3'),ylabel('2')

end
