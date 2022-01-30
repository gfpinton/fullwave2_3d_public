function [mask]=maskBdy(nX,nY,nZ,nbdy)


  mask=zeros(nX,nY,nZ);
for i=1:nX
  ri=0;
  if(i<=nbdy)
    ri=nbdy-i+1;
  end
  if(i>nX-(nbdy))
    ri=i-(nX-(nbdy));
  end

  for j=1:nY
  rj=0;
  if(j<=nbdy)
    rj=nbdy-j+1;
  end
  if(j>nY-(nbdy))
    rj=j-(nY-(nbdy));
  end
    
  for k=1:nZ
    rk=0;
    if(k<=nbdy)
      rk=nbdy-k+1;
    end
    if(k>nZ-(nbdy))
      rk=k-(nZ-(nbdy));
    end
    mask(i,j,k)=sqrt(ri^2+rj^2+rk^2);
    end
  end
end
mask=mask/max(max(max(mask)));
