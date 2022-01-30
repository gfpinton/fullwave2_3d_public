function [icvec2 icmat] = writeIcgenLagIcmat (fname,coordbin,icvec,lagmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2021-03-31
% LAST MODIFIED: 2021-04-12
% write initial conditions matrix
% if you are memory limited, use writeIcgenLag.m instead
% fname = file name, e.g. imat.dat
% coorbin = vector of size coords that writes a zero everywhere there is a zero. Put a value of one in coordbin in the active area.
% icvec = the initial condition vector as a function of time
% lagmat = lags in units of time pixels (make sure you get the sign right) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zerovec = zeros(size(icvec));

icvec2 = [zeros(1,-min(lagmat)) icvec zeros(1,max(lagmat))];
%icvecout = icvec2(lagmat(10)-min(lagmat):length(icvec)+lagmat(10)-min(lagmat));

icmat=ones(size(coordbin,1),1)*icvec2*0;

nullvec = zeros(0);
fid = fopen(fname,'wb')
fwrite(fid,nullvec);
fclose(fid)

outvec = zeros(length(coordbin),1);

fid = fopen(fname,'a')
fprintf(1,'Progress:     ');
nTic = length(icvec2);
id0 = find(coordbin<=0);
id1 = find(coordbin>0);
id1=id1(find(id1>0 & id1<=length(lagmat)));
length(id1)
max(id1)

icvec3=[icvec2 zeros(1,length(icvec2))];
for n=1:nTic
  if(length(id0))
    outvec(id0) = zerovec(n)*ones(1,length(id0));
  end
  if(length(id1))
    outvec(id1) = icvec3(lagmat(id1)-min(lagmat)+n);
  end
  fwrite(fid,outvec,'float');
  fprintf(1,'\b\b\b\b\b%0.3f',n/nTic);
  icmat(:,n)=outvec;
end
fclose(fid);
fprintf(1,'\n');



