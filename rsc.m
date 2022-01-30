  ncoordsout=size(outcoords,1);
  nRun=sizeOfFile(['genout.dat'])/4/ncoordsout; genout=readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
  figure(1), imagesc(powcompress(genout,1/4))

nX2=length(1:modX:nX); nY2=length(1:modY:nY); nZ2=length(1:modZ:nZ); 
p=reshape(genout,[],nZ2,nY2,nX2);

figure(2), imagesc(squeeze(p(end,:,round(end/2),:))), xlabel('X'), ylabel('Z'), colorbar
maxmax(genout)
minmin(genout)
sum(isnan(genout(:)))
for n=1:size(p,1)
  imagesc(squeeze(p(n,:,round(end/2),:))), xlabel('X'), ylabel('Z'), colorbar
  title(num2str(n))
  tmp=squeeze(p(n,:,:,:)); sum(isnan(tmp(:)))/length(tmp(:))
  
  drawnow
end
