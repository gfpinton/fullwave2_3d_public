clear all
addpath /nas/longleaf/home/rmjones2/path/
addpath /nas/longleaf/home/rmjones2/fullwave2_3d-master/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebecca Jones
% Updated: 12/14/21
% Matlab wrapper to launch Fullwave 2 3D code for simple neuromodulation sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1500;              % speed of sound (m/s)
omega0 = 2*pi*0.5e6;    % center radian frequency of transmitted wave
wX = 7e-2;         % lateral width of simulation field (m)
wY = 7e-2;         % elevational width of simulation field (m)
wZ = 9e-2;      % axial width of simulation field (m)
duration = wZ*2.8/c0;   % duration of simulation (s)
p0 = .6e5;              % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;               % number of points per spatial wavelength
cfl = 0.2;              % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;                % lambda
nX = round(wX/lambda*ppw);              % number of lateral elements
nY = round(wY/lambda*ppw);              % number of elevations elements
nZ=round(wZ/lambda*ppw);                % number of axial elements
nT = round(duration*c0/lambda*ppw/cfl); % number of time steps
dX = c0/omega0*2*pi/ppw;                % lateral element step size (m)
dY = c0/omega0*2*pi/ppw;                % elevational element step size(m)
dZ=dX;                                  % axial element step size (m)
dT = dX/c0*cfl;                         % time step size (s)
modX=2; modY=2; modZ=2; modT=1;         % eg. modX=2 means that your output map will contain every two lateral elements (nX of output=round(mX/modX)), modT=2 means that the simulation will record every two time steps (nT of output = round(nT/modT))
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input coordinates (transducer)

% create transducer

%focused bowl
cen=[round(nX/2) round(nY/2) (200+1)];
radh=200;
inmap=zeros(nX,nY,nZ);
for ii=1:nX
    for jj=1:nY
        for kk=1:50
            rr=(ii-cen(1))^2+(jj-cen(2))^2+((kk-(cen(3)))^2); % calculate radius
            if(rr<(radh+.5)^2 && rr>(radh-.5)^2)%&&(j-cen(1))^2/rads(1)^2+(k-cen(2))^2/rads(2)^2+(l-cen(3))^2/rads(3)^2>=0.95)
                inmap(ii,jj,kk)=1;
            end
        end
    end
end
incoords=mapToCoords3D(inmap);
figure;plot3(incoords(:,1),incoords(:,2),incoords(:,3),'.');

%create pulse
ncycles = 20; % number of cycles in pulse
dur = 10; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-2/omega0*2*pi-ncycles/2/omega0*2*pi;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
figure; plot(icvec)

icmat=ones(size(incoords,1),1)*icvec;

%for flat transducer emitting plane waves:
%{
% offset for plane
idc=1:size(incoords:,1);
for n=1:length(idc)
  i=idc(n);
  t1=t-(incoords(i,3)*dX-pcen(3)*dX)/c0;
  icvec = exp(-(1.05*t1*omega0/(ncycles*pi)).^(2*dur)).*sin(t1*omega0)*p0;
  icmat(i,:)=icvec;
end
%}

%for curved transducer or focused waves:
%calculate focal delays
fcen=cen; %focus
dd = sqrt((incoords(:,1)-fcen(1)).^2+(incoords(:,2)-fcen(2)).^2+(incoords(:,3)-fcen(3)).^2);
%calculate offset if multiple layers of transducer
idc=1:size(incoords,1);
for n=1:length(idc)
  i=idc(n);
  t1=t+(dd(i)-radh)*dX/c0;
  icvec = exp(-(1.05*t1*omega0/(ncycles*pi)).^(2*dur)).*sin(t1*omega0)*p0;
  icmat(i,:)=icvec;
end

%% Generate output coordinates (simulation field)

outcoords=coordsMatrix3d(nX,nY,nZ,modX,modY,modZ);
figure;plot3(outcoords(:,1),outcoords(:,2),outcoords(:,3),'.')
save('outcoords.mat','outcoords');

%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = zeros(nX,nY,nZ,'single')+1500; %speed of sound map %chopField(m.c,c0,orig,nX,nY);
rho = zeros(nX,nY,nZ,'single')+1000; %density map
A = zeros(nX,nY,nZ,'single'); 
beta = zeros(nX,nY,nZ,'single'); %non-linearity map
Aexp=ones(nX,nY,nZ,'single'); %attenuation map
  
% create speed of sound map with simple skull model
ROC=round(.06/dZ);%[m]
cen=[round(nX/2) round(nY/2) round(ROC*1.3)];%center of skull
skull_thickness=round(6.5e-3/dZ);
for ii=1:nX
    for jj=1:nY
        for kk=1:nZ
         rr=(ii-cen(1))^2+(jj-cen(2))^2+(kk-cen(3))^2; % calculate radius at ea>
         if(rr<(ROC+skull_thickness)^2 && rr>(ROC-0.01)^2)%&&(j-cen(1))^2/rads(>
            c(ii,jj,kk)=2800;
         end
        end
    end
end

% convert speed of sound map to density and attenuation maps
rho(c>1600)=1850;
  
A_coeff=8; % attenuation coefficient (dB//MHz/cm)
Aexp(c>1600)=dbmhzcm2aexp(A_coeff,c0,omega0,dT); 

%% launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  outdir = ['./'];
  cd(outdir)
  [Aexp2] = launch_fullwave2_3d_Aexp(c0,omega0,wX,wY,wZ,duration,p0,ppw,cfl,c,rho,A,Aexp,beta,incoords,outcoords,icmat);
  writeVabs('int',modT,'modT')
  
  ncoords=size(incoords,1);
  ncoordsout=size(outcoords,1);
  writevabs3d4('vabs3d.m',nX,nY,nZ,nT,nT,dX,dY,dZ,dT,c0,0,0,0,0,modX,modY,modZ,modT,outdir,ncoords,ncoordsout);
  
  % keep number of processes between 10 and 20
  eval('!mpirun -np 10 ./fullwave2_3d_Aexp_mpi_genout_add &') %runs in matlab window %saves output to genout.dat
  %eval('!nohup mpirun -np 10 ./fullwave2_3d_Aexp_mpi_genout_add > output2.txt &') %runs in background and saves output to output2.txt

%% process 2D slices
%run after genout.dat has finished populating (simulation has finished running)
load('outcoords')
vab3d
ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout

idc=find(outcoords(:,2)==round(nY/2)+1); %choose center slice of Y
genout_x = readGenoutSlice([outdir '/genout.dat'],1:nRun-1,size(outcoords,1),idc);

%pressure (time,Z,X)
px=reshape(genout_x,[],length(1:modZ:nZ),length(1:modX:nX));   

for ii=1:10:size(px,1)
    figure(10);imagesc(squeeze(px(ii,:,:)));
end

% intensity
pI=squeeze(sum(px).^2);
figure;imagesc((1:modX:nX),(1:modZ:nZ),dbzero(pI')',[-40 0]);colorbar;

%% find sum of pressure 3D volume

ncoordsout=size(outcoords,1)
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout; nRun=floor(nRun) 
%nRun=2000; %change nRun to load in less timeSteps

tot=zeros(ncoordsout,1)+0.0;
idc=1:ncoordsout;
genout_x = readGenoutSlice([outdir 'genout.dat'],1:nRun-1,size(outcoords,1),idc);
tot=sum(abs(genout_x(1:nRun,:,:)));

p3=reshape(tot,length(1:modZ:nZ),length(1:modY:nY),length(1:modX:nX));
figure;plot3Slice(p3); %plot 3D sum of pressure along 3 central axes
