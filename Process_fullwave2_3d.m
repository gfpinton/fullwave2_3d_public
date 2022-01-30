ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout;

idc=find(outcoords(:,2)==142);%maxcoords(1));
genout_x = readGenoutSlice([outdir 'genout.dat'],1:nRun-1,size(outcoords,1),idc>
px=reshape(genout_x,[],length(1:modZ:nZ),length(1:modX:nX));   

    %v = VideoWriter(['Pressure_propagation_c_slow.avi']);
    %open(v);
    for i=1:10:size(px,1)
        tmp=interp2easy(squeeze(px(i,:,:)),modZ,modY);
        figure(1);imagesc((1:nY)*dY*100,(1:nZ)*dZ*100,tmp)
        %imagesc(squeeze(p(i,:,:)))
        %title('Sagittal')
        ylabel('Axial Position (cm)')
        xlabel('Elevation (cm)')
        axis equal tight
        drawnow
        %frame = getframe(gcf);
        %writeVideo(v,frame);
        i
    end
    %close(v)

%pI2x=zeros(241,141,141);
%for ii=1:nRun
%pI2x=pI2x+squeeze(abs(px(ii,:,:,:)));%.^2);
%end

%px=squeeze(p(:,:,:,71));
%%%total pressure
pI2x=zeros(248,141);
for ii=1501:1531 %one cycles %1:nRun-1
pI2x=pI2x+squeeze(abs(px(ii,:,:)));%.^2);
end
%pI2x=pI2x*60e3/(1500*1000);

%%% axial pressure plot
figure;imagesc((1:nX)*dX*1000-35,(1:nZ)*dZ*1000,pI2x);axis equal tight
set(gcf,'color','w')
xlabel('Lateral Position [mm]')
ylabel('Axial Position [mm]')
axis equal tight

%%% beamplot
pvec=pI2x(8:248,71);
p_amp2=squeeze(p_amp(:,71,:)); %p_amp from Focus data
pvec_Focus=((p_amp2(:,71)));
pvec=rescale(pvec,min(pvec_Focus),max(pvec_Focus));
%pvec=interp1easy(pvec,modY);pvec=pvec-max(pvec);
%plot((1:size(pI2y,1))*dY*100,pvec)
figure;plot(pvec)
hold on;plot(pvec_Focus)
axis tight
   



