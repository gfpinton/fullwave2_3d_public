function [Aexp] = launch_fullwave2_3d_Aexp(c0,omega0,wX,wY,wZ,duration,p0,ppw,cfl,cmap,rhomap,Amap,Aexpmap,betamap,incoords,outcoords,icmat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2019-02-05
% LAST MODIFIED: 2021-03-01
% Write simulation files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optargin = size(varargin,2);

%nbdy = 20+8; % number of boundary points for PML + stencil
M=8; nbdy=40; 
lambda = c0/omega0*2*pi
				%nX = round(wX/lambda)
nX=size(cmap,1); nY=size(cmap,2); nZ=size(cmap,3);
%nXe = round(wX/lambda*ppw)+2*(nbdy+M);  % number of lateral elements
%nYe = round(wY/lambda*ppw)+2*(nbdy+M);  % number of depth elements
nXe = nX+2*(nbdy+M);  % number of lateral elements
nYe = nY+2*(nbdy+M);  % number of depth elements
nZe = nZ+2*(nbdy+M);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);

dX = c0/omega0*2*pi/ppw; dY=dX; dZ=dX;
dT = dY/c0*cfl

%Nmap=(1+boveramap/2)./(rhomap.*cmap.^4);

c=extendMap3d(cmap,nbdy+M);
rho=extendMap3d(rhomap,nbdy+M);
beta=extendMap3d(betamap,nbdy+M);
A=extendMap3d(Amap,nbdy+M);
Aexp=extendMap3d(Aexpmap,nbdy+M);
K=c.^2.*rho;
size(c);

ncoords = size(incoords,1);
ncoordsout = size(outcoords,1);
incoords = incoords+nbdy+M;
outcoords = outcoords+nbdy+M;

omega0_alpha=1e6*2*pi;

nTic=size(icmat,2);
if(optargin)
    nTic=varargin{1};
end


nu=1/2;

r=cfl;

% FOR^3D modeling 
d(2,1) = 3.26627215252963e-3 * r^7 - 7.91679373564790e-4 * r^6 + 1.08663532410570e-3 * r^5 + 2.54974226454794e-2 * r^4 + 3.23083288193913e-5 * r^3 - 3.97704676886853e-1 * r^2 + 7.95584310128586e-8 * r + 1.25425295688331;
d(3,1) = - 2.83291379048757e-3 * r^7 + 8.52796449228369e-4 * r^6 - 9.45353822586534e-4 * r^5 - 8.82015372858580e-3 * r^4 - 2.81364895458027e-5 * r^3 + 6.73021045987599e-2 * r^2 - 6.93180036837075e-8 * r - 1.23448809066664e-1;
d(4,1) = 2.32775473203342e-3 * r^7 - 5.56793042789852e-4 * r^6 + 7.77649035879584e-4 * r^5 + 2.45547234243566e-3 * r^4 + 2.31537892801923e-5 * r^3 + 1.61900960524164e-2 * r^2 + 5.70523152308121e-8 * r + 3.46683979649506e-2;
d(5,1) = - 1.68883462553539e-3 * r^7 + 3.03535823592644e-4 * r^6 - 5.64777117315819e-4 * r^5 + 2.44582905523866e-4 * r^4 - 1.68215579314751e-5 * r^3 - 2.62344345204941e-2 * r^2 - 4.14559953526389e-8 * r - 1.19918511290930e-2;
d(6,1) = 1.08994931098070e-3 * r^7 - 1.41445142143525e-4 * r^6 + 3.64794490139160e-4 * r^5 - 8.86057426195227e-4 * r^4 + 1.08681882832738e-5 * r^3 + 2.07238558666603e-2 * r^2 + 2.67876079477806e-8 * r + 4.17058420250698e-3;
d(7,1) = - 6.39950124405340e-4 * r^7 + 6.06079815415080e-5 * r^6 - 2.14633466007892e-4 * r^5 + 6.84580412267934e-4 * r^4 - 6.39907927898092e-6 * r^3 - 1.29825288653404e-2 * r^2 - 1.57775422151124e-8 * r - 1.29998325971518e-3;
d(8,1) = 2.92716539609611e-4 * r^7 - 1.87446062803024e-5 * r^6 + 9.85389372183761e-5 * r^5 - 2.40360290348543e-4 * r^4 + 2.94166215515130e-6 * r^3 + 5.57066438452790e-3 * r^2 + 7.25741366376659e-9 * r + 3.18698432679400e-4;
d(9,1) = - 6.42183857909518e-5 * r^7 + 3.38552867751042e-6 * r^6 - 2.17377151411164e-5 * r^5 + 4.98269067389945e-5 * r^4 - 6.50197868987757e-7 * r^3 - 1.19096089679178e-3 * r^2 - 1.60559948991172e-9 * r - 4.57795411807702e-5;
d(2,2) = - 4.47723278782936e-5 * r^7 - 7.69502473399932e-5 * r^6 - 1.41765498250133e-5 * r^5 - 2.54672045901272e-3 * r^4 - 4.14343385915353e-7 * r^3 + 5.00210047924752e-2 * r^2 - 1.01220354410507e-9 * r - 8.07139347787336e-8;



dmap=zeros(9,2,round(max(max(max(c))))-round(min(min(min(c)))));
for i=1:round(max(max(max(c)))-min(min(min(c))))+1
    r=((i-1)+min(min(min(c))))*dT/dX;
    
dmap(2,1,i) = 3.26627215252963e-3 * r^7 - 7.91679373564790e-4 * r^6 + 1.08663532410570e-3 * r^5 + 2.54974226454794e-2 * r^4 + 3.23083288193913e-5 * r^3 - 3.97704676886853e-1 * r^2 + 7.95584310128586e-8 * r + 1.25425295688331;
dmap(3,1,i) = - 2.83291379048757e-3 * r^7 + 8.52796449228369e-4 * r^6 - 9.45353822586534e-4 * r^5 - 8.82015372858580e-3 * r^4 - 2.81364895458027e-5 * r^3 + 6.73021045987599e-2 * r^2 - 6.93180036837075e-8 * r - 1.23448809066664e-1;
dmap(4,1,i) = 2.32775473203342e-3 * r^7 - 5.56793042789852e-4 * r^6 + 7.77649035879584e-4 * r^5 + 2.45547234243566e-3 * r^4 + 2.31537892801923e-5 * r^3 + 1.61900960524164e-2 * r^2 + 5.70523152308121e-8 * r + 3.46683979649506e-2;
dmap(5,1,i) = - 1.68883462553539e-3 * r^7 + 3.03535823592644e-4 * r^6 - 5.64777117315819e-4 * r^5 + 2.44582905523866e-4 * r^4 - 1.68215579314751e-5 * r^3 - 2.62344345204941e-2 * r^2 - 4.14559953526389e-8 * r - 1.19918511290930e-2;
dmap(6,1,i) = 1.08994931098070e-3 * r^7 - 1.41445142143525e-4 * r^6 + 3.64794490139160e-4 * r^5 - 8.86057426195227e-4 * r^4 + 1.08681882832738e-5 * r^3 + 2.07238558666603e-2 * r^2 + 2.67876079477806e-8 * r + 4.17058420250698e-3;
dmap(7,1,i) = - 6.39950124405340e-4 * r^7 + 6.06079815415080e-5 * r^6 - 2.14633466007892e-4 * r^5 + 6.84580412267934e-4 * r^4 - 6.39907927898092e-6 * r^3 - 1.29825288653404e-2 * r^2 - 1.57775422151124e-8 * r - 1.29998325971518e-3;
dmap(8,1,i) = 2.92716539609611e-4 * r^7 - 1.87446062803024e-5 * r^6 + 9.85389372183761e-5 * r^5 - 2.40360290348543e-4 * r^4 + 2.94166215515130e-6 * r^3 + 5.57066438452790e-3 * r^2 + 7.25741366376659e-9 * r + 3.18698432679400e-4;
dmap(9,1,i) = - 6.42183857909518e-5 * r^7 + 3.38552867751042e-6 * r^6 - 2.17377151411164e-5 * r^5 + 4.98269067389945e-5 * r^4 - 6.50197868987757e-7 * r^3 - 1.19096089679178e-3 * r^2 - 1.60559948991172e-9 * r - 4.57795411807702e-5;
dmap(2,2,i) = - 4.47723278782936e-5 * r^7 - 7.69502473399932e-5 * r^6 - 1.41765498250133e-5 * r^5 - 2.54672045901272e-3 * r^4 - 4.14343385915353e-7 * r^3 + 5.00210047924752e-2 * r^2 - 1.01220354410507e-9 * r - 8.07139347787336e-8;


end
r=cfl;
if(0)
%For 2D modeling:
d(2,1) = -8.74634088067635e-4 * r^7-1.80530560296097e-3 * r^6-4.40512972481673e-4 * r^5 + 4.74018847663366e-3 * r^4-1.93097802254349e-5 * r^3-2.92328221171893e-1 * r^2-6.58101498708345e-8 * r + 1.25420636437969;
d(3,1) = 7.93317828964018e-4 * r^7 + 1.61433256585486e-3 * r^6 + 3.97244786277123e-4 * r^5 + 5.46057645976549e-3 * r^4 + 1.73781972873916e-5 * r^3 + 5.88754971188371e-2 * r^2 + 5.91706982879834e-8 * r-1.23406473759703e-1;
d(4,1) = -6.50217700538851e-4 * r^7-1.16449260340413e-3 * r^6-3.24403734066325e-4 * r^5-9.11483710059994e-3 * r^4-1.41739982312600e-5 * r^3 + 2.33184077551615e-2 * r^2-4.82326094707544e-8 * r + 3.46342451534453e-2;
d(5,1) = 4.67529510541428e-4 * r^7 + 7.32736676632388e-4 * r^6 + 2.32444388955328e-4 * r^5 + 8.46419766685254e-3 * r^4 + 1.01438593426278e-5 * r^3-3.17586249260511e-2 * r^2 + 3.44988852042879e-8 * r-1.19674942518101e-2;
 d(6,1) = -2.98416281187033e-4 * r^7-3.99380750669364e-4 * r^6-1.48203388388213e-4 * r^5-6.01788793192501e-3 * r^4-6.46543538517443e-6 * r^3 + 2.41912754935119e-2 * r^2-2.19855171569984e-8 * r + 4.15554391204146e-3;
d(7,1) = 1.67882669698981e-4 * r^7 + 1.88195874702691e-4 * r^6 + 8.30579218603960e-5 * r^5 + 3.48461963201376e-3 * r^4 + 3.61873162287129e-6 * r^3-1.49875789940005e-2 * r^2 + 1.22979142197165e-8 * r-1.29213888778954e-3;
d(8,1) = -6.22209937489143e-5 * r^7-6.44890425871692e-5 * r^6-3.02936928954918e-5 * r^5-1.33386143898282e-3 * r^4-1.31215186728213e-6 * r^3 + 6.70228205200379e-3 * r^2-4.44653967516776e-9 * r + 3.15659916047599e-4;
d(9,1) = 6.84740881090240e-6 * r^7 + 1.14082245705934e-5 * r^6 + 3.03727593705750e-6 * r^5 + 2.36122782444105e-4 * r^4 + 1.26768491232397e-7 * r^3-1.53347270556276e-3 * r^2 + 4.21617557752767e-10 * r-4.51948990428065e-5;
d(2,2) = 2.13188763071246e-6 * r^7-7.41025068776257e-5 * r^6 + 2.31652037371554e-6 * r^5-2.59495924602038e-3 * r^4 + 1.20637183170338e-7 * r^3 + 5.21123771632193e-2 * r^2 + 4.42258843694177e-10 * r-4.20967682664542e-7;



dmap=zeros(9,2,round(max(max(max(c))))-round(min(min(min(c))))+1);
for i=1:round(max(max(max(c)))-min(min(min(c))))+1
    r=((i-1)+min(min(min(c))))*dT/dX;
    
    dmap(2,1,i) = -8.74634088067635e-4 * r^7-1.80530560296097e-3 * r^6-4.40512972481673e-4 * r^5 + 4.74018847663366e-3 * r^4-1.93097802254349e-5 * r^3-2.92328221171893e-1 * r^2-6.58101498708345e-8 * r + 1.25420636437969;
dmap(3,1,i) = 7.93317828964018e-4 * r^7 + 1.61433256585486e-3 * r^6 + 3.97244786277123e-4 * r^5 + 5.46057645976549e-3 * r^4 + 1.73781972873916e-5 * r^3 + 5.88754971188371e-2 * r^2 + 5.91706982879834e-8 * r-1.23406473759703e-1;
dmap(4,1,i) = -6.50217700538851e-4 * r^7-1.16449260340413e-3 * r^6-3.24403734066325e-4 * r^5-9.11483710059994e-3 * r^4-1.41739982312600e-5 * r^3 + 2.33184077551615e-2 * r^2-4.82326094707544e-8 * r + 3.46342451534453e-2;
dmap(5,1,i) = 4.67529510541428e-4 * r^7 + 7.32736676632388e-4 * r^6 + 2.32444388955328e-4 * r^5 + 8.46419766685254e-3 * r^4 + 1.01438593426278e-5 * r^3-3.17586249260511e-2 * r^2 + 3.44988852042879e-8 * r-1.19674942518101e-2;
 dmap(6,1,i) = -2.98416281187033e-4 * r^7-3.99380750669364e-4 * r^6-1.48203388388213e-4 * r^5-6.01788793192501e-3 * r^4-6.46543538517443e-6 * r^3 + 2.41912754935119e-2 * r^2-2.19855171569984e-8 * r + 4.15554391204146e-3;
dmap(7,1,i) = 1.67882669698981e-4 * r^7 + 1.88195874702691e-4 * r^6 + 8.30579218603960e-5 * r^5 + 3.48461963201376e-3 * r^4 + 3.61873162287129e-6 * r^3-1.49875789940005e-2 * r^2 + 1.22979142197165e-8 * r-1.29213888778954e-3;
dmap(8,1,i) = -6.22209937489143e-5 * r^7-6.44890425871692e-5 * r^6-3.02936928954918e-5 * r^5-1.33386143898282e-3 * r^4-1.31215186728213e-6 * r^3 + 6.70228205200379e-3 * r^2-4.44653967516776e-9 * r + 3.15659916047599e-4;
dmap(9,1,i) = 6.84740881090240e-6 * r^7 + 1.14082245705934e-5 * r^6 + 3.03727593705750e-6 * r^5 + 2.36122782444105e-4 * r^4 + 1.26768491232397e-7 * r^3-1.53347270556276e-3 * r^2 + 4.21617557752767e-10 * r-4.51948990428065e-5;
dmap(2,2,i) = 2.13188763071246e-6 * r^7-7.41025068776257e-5 * r^6 + 2.31652037371554e-6 * r^5-2.59495924602038e-3 * r^4 + 1.20637183170338e-7 * r^3 + 5.21123771632193e-2 * r^2 + 4.42258843694177e-10 * r-4.20967682664542e-7;
end
r=cfl;
end


dcmap=round(c)-min(min(min(c)))+1;
dcmap=round(c)-min(min(min(c)));

%% PMLS %%
kpml=1;
dpmlxOld=zeros(1,nXe); alphapmlxOld=dpmlxOld; apmlxOld=dpmlxOld; bplm=dpmlxOld; 
dpmlxOld=zeros(1,nXe); alphapmlxOld=dpmlxOld; apmlxOld=dpmlxOld; bpmlxOld=dpmlxOld+1;
dpmlxOld2=zeros(1,nXe); alphapmlxOld2=dpmlxOld; apmlxOld2=dpmlxOld; bpmlxOld2=dpmlxOld+1;
L=dX*nbdy;
Rc=1e-30;
d0=-3*c0*log(Rc)/(2*L)
for i=1:nbdy
    dpmlxOld(i+(nXe-M-nbdy+1))=d0*(i/nbdy)^2;
    dpmlxOld(M+nbdy+1-i)=d0*(i/nbdy)^2;
end
for i=1:nbdy
    alphapmlxOld(i+(nXe-M-nbdy+1))=omega0/2*(nbdy-i)/nbdy;
    alphapmlxOld(M+nbdy+1-i)=omega0/2*(nbdy-i)/nbdy;
end
alphapmlxOld=alphapmlxOld/10;
%bpmlxOld=exp(-(dpmlxOld/kpml+alphapmlxOld)*dT);
%apmlxOld=dpmlxOld/(kpml*(dpmlxOld+kpml*alphapmlxOld))*(bpmlxOld-1);
[apmlxOld bpmlxOld]=ab(dpmlxOld,kpml,alphapmlxOld,dT);

for i=1:nbdy
    dpmlxOld2(i+(nXe-M-nbdy+1))=d0*((i-1/2)/nbdy)^2;
    dpmlxOld2(M+nbdy+1-i-1)=d0*((i-1/2)/nbdy)^2;
end
for i=1:nbdy
    alphapmlxOld2(i+(nXe-M-nbdy+1))=omega0/2*(nbdy-(i-1/2))/nbdy;
    alphapmlxOld2(M+nbdy+1-i-1)=omega0/2*(nbdy-(i-1/2))/nbdy;
end
alphapmlxOld2=alphapmlxOld2/10;
%bpmlxOld2=exp(-(dpmlxOld2/kpml+alphapmlxOld2)*dT);
%apmlxOld2=dpmlxOld2/(kpml*(dpmlxOld2+kpml*alphapmlxOld2))*(bpmlxOld2-1);
[apmlxOld2 bpmlxOld2]=ab(dpmlxOld2,kpml,alphapmlxOld2,dT);

plot(apmlxOld), hold on
plot(apmlxOld2,'r')
%plot(alphapmlxOld2,'g')
hold off

plot(bpmlxOld), hold on
plot(bpmlxOld2,'r')
hold off

%%

dpmlyOld=zeros(1,nYe); alphapmlyOld=dpmlyOld; apmlyOld=dpmlyOld; bplm=dpmlyOld; 
dpmlyOld=zeros(1,nYe); alphapmlyOld=dpmlyOld; apmlyOld=dpmlyOld; bpmlyOld=dpmlyOld+1; 
dpmlyOld2=zeros(1,nYe); alphapmlyOld2=dpmlyOld; apmlyOld2=dpmlyOld; bpmlyOld2=dpmlyOld+1;

for i=1:nbdy
    dpmlyOld(i+(nYe-M-nbdy+1))=d0*(i/nbdy)^2;
    dpmlyOld(M+nbdy+1-i)=d0*(i/nbdy)^2;
end
for i=1:nbdy
    alphapmlyOld(i+(nYe-M-nbdy+1))=omega0/2*(nbdy-i)/nbdy;
    alphapmlyOld(M+nbdy+1-i)=omega0/2*(nbdy-i)/nbdy;
end
alphapmlyOld=alphapmlyOld/10;
%bpmlyOld=exp(-(dpmlyOld/kpml+alphapmlyOld)*dT);
%apmlyOld=dpmlyOld/(kpml*(dpmlyOld+kpml*alphapmlyOld))*(bpmlyOld-1);
[apmlyOld bpmlyOld]=ab(dpmlyOld,kpml,alphapmlyOld,dT);

for i=1:nbdy
    dpmlyOld2(i+(nYe-M-nbdy+1))=d0*((i-1/2)/nbdy)^2;
    dpmlyOld2(M+nbdy+1-i-1)=d0*((i-1/2)/nbdy)^2;
end
for i=1:nbdy
    alphapmlyOld2(i+(nYe-M-nbdy+1))=omega0/2*(nbdy-(i-1/2))/nbdy;
    alphapmlyOld2(M+nbdy+1-i-1)=omega0/2*(nbdy-(i-1/2))/nbdy;
end
alphapmlyOld2=alphapmlyOld2/10;
%bpmlyOld2=exp(-(dpmlyOld2/kpml+alphapmlyOld2)*dT);
%apmlyOld2=dpmlyOld2/(kpml*(dpmlyOld2+kpml*alphapmlyOld2))*(bpmlyOld2-1);
[apmlyOld2 bpmlyOld2]=ab(dpmlyOld2,kpml,alphapmlyOld2,dT);

%apmlyOld=apml; bpmlyOld=bpml; apmlyOld2=apml2; bpmlyOld2=bpml2;


plot(apmlyOld), hold on
plot(apmlyOld2,'r')
hold off

plot(bpmlyOld), hold on
plot(bpmlyOld2,'r')
hold off

%%

dpmlzOld=zeros(1,nZe); alphapmlzOld=dpmlzOld; apmlzOld=dpmlzOld; bplm=dpmlzOld; 
dpmlzOld=zeros(1,nZe); alphapmlzOld=dpmlzOld; apmlzOld=dpmlzOld; bpmlzOld=dpmlzOld+1; 
dpmlzOld2=zeros(1,nZe); alphapmlzOld2=dpmlzOld; apmlzOld2=dpmlzOld; bpmlzOld2=dpmlzOld+1;

for i=1:nbdy
    dpmlzOld(i+(nZe-M-nbdy+1))=d0*(i/nbdy)^2;
    dpmlzOld(M+nbdy+1-i)=d0*(i/nbdy)^2;
end
for i=1:nbdy
    alphapmlzOld(i+(nZe-M-nbdy+1))=omega0/2*(nbdy-i)/nbdy;
    alphapmlzOld(M+nbdy+1-i)=omega0/2*(nbdy-i)/nbdy;
end
alphapmlzOld=alphapmlzOld/10;
%bpmlzOld=exp(-(dpmlzOld/kpml+alphapmlzOld)*dT);
%apmlzOld=dpmlzOld/(kpml*(dpmlzOld+kpml*alphapmlzOld))*(bpmlzOld-1);
[apmlzOld bpmlzOld]=ab(dpmlzOld,kpml,alphapmlzOld,dT);

for i=1:nbdy
    dpmlzOld2(i+(nZe-M-nbdy+1))=d0*((i-1/2)/nbdy)^2;
    dpmlzOld2(M+nbdy+1-i-1)=d0*((i-1/2)/nbdy)^2;
end
for i=1:nbdy
    alphapmlzOld2(i+(nZe-M-nbdy+1))=omega0/2*(nbdy-(i-1/2))/nbdy;
    alphapmlzOld2(M+nbdy+1-i-1)=omega0/2*(nbdy-(i-1/2))/nbdy;
end
alphapmlzOld2=alphapmlzOld2/10;
%bpmlzOld2=exp(-(dpmlzOld2/kpml+alphapmlzOld2)*dT);
%apmlzOld2=dpmlzOld2/(kpml*(dpmlzOld2+kpml*alphapmlzOld2))*(bpmlzOld2-1);
[apmlzOld2 bpmlzOld2]=ab(dpmlzOld2,kpml,alphapmlzOld2,dT);

%apmlzOld=apml; bpmlzOld=bpml; apmlzOld2=apml2; bpmlzOld2=bpml2;


plot(apmlzOld), hold on
plot(apmlzOld2,'r')
hold off

plot(bpmlzOld), hold on
plot(bpmlzOld2,'r')
hold off


Amask=1-maskBdy(size(Aexp,1),size(Aexp,2),size(Aexp,3),nbdy+M);
%Amask=Amask.^2;
Aexp=Aexp.*Amask;

%%

writeMapXYZ('c.dat',c);
writeMapXYZ('K.dat',K);
writeMapXYZ('rho.dat',rho);
writeMapXYZ('beta.dat',beta);
writeMapXYZ('Aexp.dat',Aexp);

writeCoords('icc.dat',incoords-1);
writeCoords('outc.dat',outcoords-1);
writeIC('icmat.dat',icmat');

writeVabs('float',dX,'dX',dY,'dY',dZ,'dZ',dT,'dT',c0,'c0');
writeVabs('int',nXe,'nX',nYe,'nY',nZe,'nZ',nT,'nT',ncoords,'ncoords',ncoordsout,'ncoordsout',nTic,'nTic');


fid = fopen(['d.dat'],'wb'); fwrite(fid,d','float'); fclose(fid);

fid = fopen(['dmap.dat'],'wb'); fwrite(fid,zeros(0),'float'); fclose(fid);

fid = fopen(['dmap.dat'],'ab'); 
fwrite(fid,zeros(0),'float'); 
for i=1:size(dmap,1)
    for j=1:size(dmap,2)
        for k=1:size(dmap,3)
            fwrite(fid,dmap(i,j,k),'float');
        end
    end
end
fclose(fid);

ndmap=size(dmap,3)
writeVabs('int',ndmap,'ndmap');

writeMapXYZ('dcmap.dat',dcmap,'int');

