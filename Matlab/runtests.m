% run tests for the stacking functions
close all;
load('testdata/Xcorr_remEQ_xcorr_UW.HOOD.ZH.ALS0.ZZ.mat');

%% plot raw data
[nsamp,ntrace]=size(xcorrdata);
tt=timeflag(:,1);
data_norm=xcorrdata;
for i=1:ntrace
    data_norm(:,i)=xcorrdata(:,i)./max(abs(xcorrdata(:,i)));
end
% hold off;
%
data_norm2=data_norm';


dt=abs(tt(2)-tt(1));
zeroerror=dt/100; % used to define whether the tt is zero, bounded by -1*thisvalue and 1*thisvalue.
zeroidx=find(tt > -1.0*zeroerror & tt < zeroerror);

vmin=.7; vmax=4.5;
taxisext=1.5; %extension of time axis. time axis in plotting will be taxisext*tmax0;
tmin0=metadata.DIST/vmax;
tmax0=100+metadata.DIST/vmin;
if tmax0 > tt(end)
    tmax0 = tt(end);
end

tmin=-taxisext*tmax0;
tmax=taxisext*tmax0;
if tmin < tt(1); tmin=tt(1);end
if tmax > tt(end); tmax=tt(end);end
itmin=round((tmin - tt(1))/dt)+1;
itmax=round((tmax - tt(1))/dt)+1;

noiseoffset= 50; %noise window will start after signal window with this 'noiseoffset' in seconds.
noiseoffsetnpts=round(noiseoffset/dt);
clear idx0 idx00;
idx0=find(tt > tmin0-.5*dt & tt < tmin0+.5*dt);
idx00=find(tt > tmax0-.5*dt & tt < tmax0+.5*dt);
signalwin_pos=[idx0,idx00];
noisewin_pos=signalwin_pos + noiseoffsetnpts + idx00 - idx0 + 1;
if noisewin_pos(2)>nsamp
    noisewin_pos(2)=nsamp;
end
%negative 
signalwin_neg=[2*zeroidx-signalwin_pos(2),2*zeroidx-signalwin_pos(1)];
noisewin_neg=[2*zeroidx-noisewin_pos(2),2*zeroidx-noisewin_pos(1)];
%% stacking
nstack=9;
par0=struct('verbose',0);
timestackall=nan(nstack,1);
%linear
[dmean,statmean]=seisstack(xcorrdata,'linear');
timestackall(1)=statmean.t;

%linear - normalized before stacking
lpar.normalize=1;
[dmean_norm,statmean_norm]=seisstack(xcorrdata,'linear',lpar);
timestackall(2)=statmean_norm.t;

%selective - ccmin=0.0
[dselect,statselective]=seisstack(xcorrdata,'selective');
timestackall(3)=statselective.t;

%robust
% rpar.reference=pwstack(xcorrdata);
[drobust,statrobust]=seisstack(xcorrdata,'robust');
timestackall(4)=statrobust.t;

%robust - specified window for weight calculation
robustitmin=signalwin_neg(1);
robustitmax=signalwin_pos(2);
rpar.win=[robustitmin,robustitmax];
[drobust_win,statrobust_win]=seisstack(xcorrdata,'robust',rpar);
timestackall(5)=statrobust_win.t;  

%pws
[dpws,statpws]=seisstack(xcorrdata,'pws');
timestackall(6)=statpws.t;

%tf-pws
par=par0;
par.pow=0.5;
[dtfpws,stattfpws]=seisstack(xcorrdata,'tf-pws',par);
timestackall(7)=stattfpws.t;

%N^th root stacking
par=par0;
par.N=2;
[dnroot,statnroot]=seisstack(xcorrdata,'nroot',par);
timestackall(8)=statnroot.t;

%acfstack
par=par0;
par.window=round(200/dt);
par.overlap=0.9;
par.harshness=0.2;
[dacf,statacf]=seisstack(xcorrdata,'acf',par);
timestackall(9)=statacf.t;

dstackall=[dmean, dmean_norm, dselect,drobust,drobust_win,...
    dpws,dtfpws,dnroot,dacf];%
%%
figure('Position',[400 400 1200 600]); 
subplot(1,2,1);hold on;
hi=image(tt,1:ntrace,data_norm2);
hi.CDataMapping='scaled';
pdatamask=ones(size(data_norm2));
pdatamask(isnan(data_norm2))=0;
hi.AlphaData=pdatamask;

colormap('jet');
plot([0 0],[1 ntrace],'k','LineWidth',1)
plot([tt(signalwin_pos(1)) tt(signalwin_pos(1))],[1 ntrace],'k--','linewidth',1.5);
plot([tt(signalwin_pos(2)) tt(signalwin_pos(2))],[1 ntrace],'k--','linewidth',1.5);

plot([tt(signalwin_neg(1)) tt(signalwin_neg(1))],[1 ntrace],'k--','linewidth',1.5);
plot([tt(signalwin_neg(2)) tt(signalwin_neg(2))],[1 ntrace],'k--','linewidth',1.5);

text(tt(round(mean(signalwin_pos)))+10,0.9*ntrace,'signal','HorizontalAlignment','center','fontsize',12);
text(tt(round(mean(signalwin_pos)))+10,0.84*ntrace,'window','HorizontalAlignment','center','fontsize',12);

text(tt(round(mean(signalwin_neg))),0.9*ntrace,'signal','HorizontalAlignment','center','fontsize',12);
text(tt(round(mean(signalwin_neg))),0.84*ntrace,'window','HorizontalAlignment','center','fontsize',12);
hold off;
hc=colorbar;
hc.Label.String='Normalized amplitude';
hc.Location='eastoutside';
set(hc,'FontSize',12,'TickDir','out');
caxis([-1 1]);
% xlim([tt(1) tt(end)])
xlim([tmin tmax])
ylim([1 ntrace])
title('Raw NCFs');

xlabel('Cross-correlation time (s)')
axis on;
box on;
set(gca,'TickDir','out','YDir','norm','Fontsize',12)
drawnow;

% plot stacks
stackscale=0.4;
lwidth=1;
% cmap=cool(2*nstack);
cmap=[0.0 0.0 0.0;
      0.4 0.4 0.4;
      0.0 0.5 0.6;
      0.0 0.8 0.8;
      0.0 0.0 1.0;
      0.0 0.5 0.0;
      0.0 0.9 0.0;
      1.0 0.0 0.0;
      1.0 0.0 1.0];
% figure('Position',[400 400 700 700]);
subplot(1,2,2);
hold on;
for i=1:nstack
    plot(tt,i+stackscale*dstackall(:,i)/max(abs(dstackall(itmin:itmax,i))),'k','linewidth',lwidth,'color',cmap(i,:));
end

xlim([tmin tmax]);
ylim([0.5 nstack+0.5]);
plot([0 0],[0.5 nstack+0.5],'k');

%text label
text(-0.95*taxisext*tmax0,.7,'Linear stack of all traces','fontsize',12);
text(-0.95*taxisext*tmax0,1.7,'Linear stack after normalization','fontsize',12);
text(-0.95*taxisext*tmax0,2.7,['Selective with threshold: ' num2str(statselective.par.ccmin)],'fontsize',12);
text(-0.95*taxisext*tmax0,3.7,'Robust window: entire trace','fontsize',12);
text(-0.95*taxisext*tmax0,4.7,['Robust window: ',num2str(tt(robustitmin)),' s to ',num2str(tt(robustitmax)),' s'],'fontsize',12);
text(-0.95*taxisext*tmax0,5.7,['Phase-weighted stack (power=',num2str(statpws.par.pow),')'],'fontsize',12);
text(-0.95*taxisext*tmax0,6.7,['Time-frequency PWS (power=',num2str(stattfpws.par.pow),')'],'fontsize',12);
text(-0.95*taxisext*tmax0,7.7,['N^t^h root (N=',num2str(statnroot.par.N),')'],'fontsize',12);
text(-0.95*taxisext*tmax0,8.7,'Adaptive covariance filter','fontsize',12);

for i=1:nstack
    text(0.9*tmax0,i-0.3,['t=',num2str(timestackall(i)), ' s'],'fontsize',12);
end

hold off;
axis on;
box on;
grid on;
title('Comparison of stacks');
set(gca,'FontSize',14,'YDir','reverse','YTick',1:nstack,'YTickLabel',...
    {'Linear-1','Linear-2','Selective','Robust-1','Robust-2','PWS','tf-PWS','N^t^h root','ACF'});
drawnow;

set(gcf,'PaperPositionMode','auto');   
eval(['print -dpng -r300 comparison.png']);

%%
zeroidx=find(tt>=-0.001*dt & tt<=0.001*dt);
fs=1/dt;
freqinc=fs/400;
dstackallnorm_neg=nan(zeroidx,nstack);
dstackallnorm_pos=nan(zeroidx,nstack);
for i=1:nstack
    dstackallnorm_neg(:,i)=dstackall(1:zeroidx,i)/max(abs(dstackall(1:zeroidx,i)));
    dstackallnorm_pos(:,i)=dstackall(zeroidx:end,i)/max(abs(dstackall(zeroidx:end,i)));
end
% negative side
figure('Position',[400 400 1000 450]);
[px,fx]=pspectrum(dstackallnorm_neg,fs,'FrequencyResolution',freqinc);
pdb=pow2db(px);
pdb(isinf(pdb))=nan;

subplot(1,2,1);hold on;
for i=1:nstack
    lw=lwidth;
    if i==1;lw=2*lw;end
    plot(fx,pdb(:,i),'k','linewidth',lw,'color',cmap(i,:));
end
xlim([0.01 max(fx)]);
ylim([-160 0]);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Power spectrum (dB)','FontSize',14);
legend('Linear-1','Linear-2','Selective','Robust-1','Robust-2','PWS','tf-PWS','N^t^h root','ACF',...
    'Location','south');
title('(a) Negative lags');
set(gca,'fontsize',14,'XScale','log');
hold off;
axis on; box on;grid on;

[px,fx]=pspectrum(dstackallnorm_pos,fs,'FrequencyResolution',freqinc);
pdb=pow2db(px);
pdb(isinf(pdb))=nan;

subplot(1,2,2);hold on;
for i=1:nstack
    lw=lwidth;
    if i==1;lw=2*lw;end
    plot(fx,pdb(:,i),'k','linewidth',lw,'color',cmap(i,:));
end
xlim([.01 max(fx)]);
ylim([-160 0]);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Power spectrum (dB)','FontSize',14);
legend('Linear-1','Linear-2','Selective','Robust-1','Robust-2','PWS','tf-PWS','N^t^h root','ACF',...
    'Location','south');
title('(b) Positive lags');
set(gca,'fontsize',14,'XScale','log');
hold off;
axis on; box on; grid on;

set(gcf,'PaperPositionMode','auto');   
eval(['print -dpng -r300 comparison_spectrum.png']);