function [] = writeMap(fname,cmap,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2012-01-23
% LAST MODIFIED: 2013-11-13
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
for j=1:size(cmap,2)
  fprintf(1,'\b\b\b\b\b%0.3f',j/size(cmap,2));
  for i=1:size(cmap,1)
    vec = cmap(i,j,:); fwrite(fidc,vec,strtyp);
  end
end
fclose(fidc);
fprintf(1,'\n');
