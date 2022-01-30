function [] = writeMapXYZ(fname,cmap,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2012-01-23
% LAST MODIFIED: 2021-03-01
% write map to disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  strtyp='float';
   optargin = size(varargin,2);
  if(optargin==1)
    strtyp=varargin{1}
  end
fprintf(1,'Progress:     ');
fid = fopen(fname,'wb'); fwrite(fid,zeros(0)); fclose(fid);
fidc = fopen(fname,'a');
for i=1:size(cmap,1)
  fprintf(1,'\b\b\b\b\b%0.3f',i/size(cmap,1));
  for j=1:size(cmap,2)
    vec = cmap(i,j,:); fwrite(fidc,vec,strtyp);
  end
end
fclose(fidc);
fprintf(1,'\n');
