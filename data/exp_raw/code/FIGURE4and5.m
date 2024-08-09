clear all;clc;close all;
load Px.mat;
conName={'normal sleep','sleep deprivation'};

%% Step 1:
PxCz=reshape(Px(:,13,:),[4097,2,66]);
Pxm=mean(PxCz,3);
Pxs=std(PxCz,0,3)./sqrt(66);
xfreq=[0:4096]./4096*250;
figure(1)

ERP=Pxm;

ERPse=Pxs;
xtime=[0:4096]./4096*250;
colorM=[1 0 1;0 1 0;];
colorS=[1 0.8 1;0.8 1 0.8];

ERP1=ERP+ERPse;
ERP2=ERP-ERPse;
ERP12=[ERP1;ERP2(size(ERP2,1):-1:1,:)];
xtime12=[xtime xtime(end:-1:1)];
for i=1:2
    patch(xtime12,ERP12(:,i),colorS(i,:),'EdgeColor','none');
    hold on;
end
for i=1:2
    plot(xtime,ERP(:,i),'color',colorM(i,:),'LineWidth',3);
    hold on;
end
xlim([2 40]);
ylim([-32 10])
grid on;
legend(conName)

%% Step 2:
%%%% theta (4–8 Hz), alpha (8–13 Hz), beta (13–30 Hz)
freqRang={[67:132],[133:213],[214:492]};
for i=1:3
    Rhythm(:,:,i)=squeeze(mean(Px(freqRang{i},:,:),1));%  61 66*2 3
end
Rhythm=reshape(Rhythm,[61 2 66 3]);
Rhythmm=squeeze(mean(Rhythm,3));% 61 2 3

load chanlocs.mat;
rhthName={'theta','alpha','beta'};
figure(2)
for i=1:2
    for j=1:3
        subplot(2,3,i*3+j-3)
        topoplot(squeeze(Rhythmm(:,i,j)),chanlocs,'maplimits',[-20 10]);
        title([conName{i},'-',rhthName{j}]);
        colorbar;
    end
end



