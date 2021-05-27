clear
load('SingleExcitPYRPVData.mat')
Single = META;
clearvars -EXCEPT Single 

load('PairedExcitPYRPVData.mat')
Paired = GOODPairedMETA;
ALLPairedMETA_PV = ALLPairedMETA;
clearvars -EXCEPT Single Paired ALLPairedMETA_PV;
GoodGreen = [0 0.7 0];
GoodCyan = [0 0.9 0.9];


samplePV     = 10; %My beautiful PV <3
sampleWT     = 9; %Didn't really like this example but Dean liked it so...
AdaptCurrent = 8;

%% Getting Data from Individually Recorded Cells
% Getting Indexes of ProperVariables
T1 = find([Single.Group]==1); %PV Indexes
T2 = find([Single.Group]==2); %WT Indexes
% TargetInt1 = find(Single(T1(1)).CONDINT==[NaN NaN 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4]);
% TargetInt2 = find(Single(T2(1)).CONDINT==[NaN NaN 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4]);

%Calculating Variables for Spike Num
SingleSpikesPV = cell2mat({Single(T1).Spikenum}); SingleSpikesPV = SingleSpikesPV(2:10,:)';
SingleSpikesPYR = cell2mat({Single(T2).Spikenum}); SingleSpikesPYR = SingleSpikesPYR(2:10,:)';

% %Calculating ISI for the example trace
% ISIWT = diff(Paired(sampleWT).SpData{1,8}(:,2)'); %from the sample trace at intensity 0.3nA
% ISIPV = diff(Single(samplePV).SpData{8}(:,2)'); %from the sample trace at intensity 0.3nA

%Calculating Average Adaptation Index
AdaptIndexPV = cell2mat({Single(T1).AdaptIndex}); AdaptIndexPV = AdaptIndexPV(8,:)';
AdaptIndexPYR = cell2mat({Single(T2).AdaptIndex}); AdaptIndexPYR = AdaptIndexPYR(8,:)';
SingleAdaptIndex = {AdaptIndexPYR,AdaptIndexPV};

%% Getting Data from Simultaneous Recording
BadVarName    = squeeze(struct2cell(Paired)); 
Ind10Steps    = cellfun(@length,BadVarName(3,:))>9; %Finding Files with 10 or more Excitability current steps
Ind13Steps    = cellfun(@length,BadVarName(3,:))>12; %Finding Files with 13 or more Excitability current steps
Paired10Steps = Paired(Ind10Steps);
Paired13Steps = Paired(Ind13Steps);

PairedAdaptIndex      = zeros(length(Paired10Steps),2);
PairedSpikenumPYR     = zeros(length(Paired10Steps),9);
PairedSpikenumPV      = zeros(length(Paired13Steps),12);

for d = 1:length(Paired10Steps)
    PairedAdaptIndex(d,:)  = Paired10Steps(d).AdaptIndex(AdaptCurrent,:);
    PairedSpikenumPYR(d,:) = Paired10Steps(d).SpikenumPYR(2:10)'; %We don't count the first one because it is a negative current
end

for d = 1:length(Paired13Steps)
    PairedSpikenumPV(d,:)  = Paired13Steps(d).SpikenumPV(2:13)';
end

%% Combining Single and Simultaneous Data
AdaptIndex = {[SingleAdaptIndex{1};PairedAdaptIndex(:,1)],[SingleAdaptIndex{2};PairedAdaptIndex(:,2)]};

SpikesPV = PairedSpikenumPV;
SpikesPYR = [SingleSpikesPYR;PairedSpikenumPYR];

AvrSpikesPYR = mean(SpikesPYR);
SEMPYR = std(SpikesPYR)/sqrt(size(SpikesPYR,1));

%The N for the last 4 values in the SEMPV is different than the rest. 
AvrSpikesPV = mean(SpikesPV);
SEMPV =  std(SpikesPV)/sqrt(size(SpikesPV,1));

XPYR = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];
XPV = [XPYR 0.5 0.55 0.6];
%% Plotting Figure 1
close all

MainFigure = figure;
set(gcf,'color','w','position',[80  60  1160  1680]);

SP1 = subplot('Position',[0.03 0.77 0.28 0.23]);
SP2 = subplot('Position',[0.36 0.77 0.28 0.23]);
SP3a = subplot('Position',[0.06 0.4 0.21 0.3 ]);
SP3b = subplot('Position',[0.28 0.43 0.02 0.27]);
SP4a = subplot('Position',[0.37 0.4 0.21 0.3]);
SP4b = subplot('Position',[0.60 0.43 0.02 0.27]);
% SP5 = subplot('Position',[0.04 0.05 0.24 0.29]);
% SP6 = subplot('Position',[0.38 0.05 0.24 0.29]);
AX1 = axes('Position',[0 0 1 1],'Visible','off');

SP7 = subplot('Position',[0.68 0.77 0.28 0.23]);
SP8a = subplot('Position',[0.705 0.4 0.21 0.3]);
SP8b = subplot('Position',[0.9228 0.43 0.065 0.27]);
% SP9 = subplot('Position',[0.73 0.05 0.24 0.29]);
% AX2 = axes('Position',[0.09 0.43 0.296 0.3]);
% AX3 = axes('Position',[0.55 0.41 0.296 0.3]);

%%
% Plotting Raw Traces
subplot(SP1)
set(gcf,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; .75 0 .75; .75 .75 0; .25 .25 .25])
x1 = (0:3500)';
zWT = (Paired(sampleWT).Data{1}([1 3 5 9],1500:5000)');
xMat = repmat(x1,1,size(zWT,2));
zMatWT = zWT;
plot(xMat,zMatWT,'linewidth',3)
ylim([-80 20])
% legend('-0.1 nA','0.1 nA','0.2 nA','0.4 nA','FontSize',12,'FontWeight','bold')
% legend('boxoff')
box off
axis off

%Crazy 3D stuff
% yMat = repmat((1:10),length(x1),1);
% plot3(xMat,yMat,zMat,'linewidth',2.5)
% view(40,40)
%print -djpeg100 PVTraces3D.jpg -r300

subplot(SP2)
set(gcf,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; .75 0 .75; .75 .75 0; .25 .25 .25])
zPV = (Single(samplePV).Data([1 3 5 9],1500:5000)');
zMatPV = zPV;
plot(xMat,zMatPV,'linewidth',3)
ylim([-75 30])
box off
axis off

% %Plotting ISI and Adaptation Index
% subplot(SP3a)
% plot(1:length(ISIWT),ISIWT,'Color',GoodGreen,'Marker','o','LineWidth',4,'MarkerFaceColor',GoodGreen)
% box off
% set(gca,'linewidth',4,'fontsize',16,'fontweight','bold','XTick',[1 2 3 4 5])
% title('Inter-Spike Interval Pyr','FontSize',18)
% SP3a.Title.Position = [3.5 73.7883 0];
% xlabel('ISI Number','FontSize',18)
% ylabel('ISI(ms)','FontSize',18)
% xlim([0.75 length(ISIWT)+0.25])
% ylim([0 70])
% 
% subplot(SP3b)
% boxplot(AdaptIndex{1},'PlotStyle','compact','Colors',GoodGreen,'Labels',{''});
% h = gca;
% set(h,'linewidth',4,'fontsize',16,'fontweight','bold','yaxislocation','right','yticklabel',[])
% h.Children.Children(6).LineWidth = 2;
% h.Children.Children(5).LineWidth = 8;
% h.XAxis.Visible = 'off';
% ylim([-0.025 0.16])
% xlim([0.75 1.25])
% box off
% 
% subplot(SP4a)
% plot(1:length(ISIPV),ISIPV,'Color','r','Marker','o','LineWidth',4,'MarkerFaceColor','r')
% box off
% set(gca,'linewidth',4,'fontsize',16,'fontweight','bold')
% title('Inter-Spike Interval PV','FontSize',18)
% SP4a.Title.Position = [9.5 73.7883 0];
% xlabel('ISI Number','FontSize',18)
% xlim([0.5 length(ISIPV)+0.25])
% ylim([0 70])
% 
% subplot(SP4b)
% boxplot(AdaptIndex{2},'PlotStyle','compact','Colors','r','Labels',{''})
% h = gca;
% set(h.Children.Children(2),'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
% set(h,'linewidth',4,'fontsize',16,'fontweight','bold','yaxislocation','right','yticklabel',[])
% h.Children.Children(6).LineWidth = 2;
% h.Children.Children(5).LineWidth = 8;
% h.XAxis.Visible = 'off';
% % ylabel('Adaptation Idex','FontSize',12)
% ylim([-0.025 0.16])
% box off

% % Plotting SpikeNums
% subplot(SP5)
% errorbar(XPYR,AvrSpikesPYR,SEMPYR,'color','k','LineWidth',4,'Marker','o','MarkerSize',4)
% set(gca,'fontsize',10,'fontweight','bold')
% title('Input/Output PYR','FontSize',15)
% xlabel('Input Intensity (nA)','FontSize',12)
% ylabel('Spk Count','FontSize',12)
% xlim([XPYR(1)*0.95 XPYR(end)*1.05])
% ylim([0 18])
% box off
% 
% subplot(SP6)
% errorbar(XPV,AvrSpikesPV,SEMPV,'color','r','LineWidth',4,'Marker','o','MarkerSize',4)
% hold on
% % scatter([0.5 0.55 0.6],[14 18 20],55,[1 0 0],'d','MarkerFaceColor','r')
% hold off
% box off
% set(gca,'fontsize',10,'fontweight','bold')
% title('Input/Output PV','FontSize',15)
% xlabel('Input Intensity (nA)','FontSize',12)
% ylabel('Spk Count','FontSize',12)
% xlim([XPV(1)*0.95 XPV(end)*1.05])
% ylim([0 18])
% box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SST Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear AdaptIndex BadVarName Ind10Steps Ind13Steps  
load('PairedExcitSSTData.mat')
SSTMETA = GOODPairedMETA;
ALLPairedMETA_SST = ALLPairedMETA;

sampleSST = 14; %I think people will like this one
ISISST = diff(SSTMETA(sampleSST).spData{9,2}); %from the sample trace at intensity 0.3nA

BadVarName = squeeze(struct2cell(SSTMETA)); 
Ind10Steps = find(cellfun(@length,BadVarName(3,:))>9); %Finding Files with 10 or more Excitability current steps
SST10Steps = SSTMETA(Ind10Steps);

for j = 1:length(SST10Steps)
    AdaptIndex(j,:)   = SST10Steps(j).AdaptIndex(AdaptCurrent,2);
    SpikesSST(j,:)     = SST10Steps(j).SpikenumSST(2:10);
end

% for j = 1:length(SST13Steps) 
%     SpikesSSTLater(j,:) = SST13Steps(j).SpikenumSST(11:13);
% end

AvrSpikesSST = mean(SpikesSST);
SEMSST       =  std(SpikesSST)/sqrt(size(SpikesSST,1));
XSST         = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];

%% Plotting SST

subplot(SP7)
set(gcf,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; .75 0 .75; .75 .75 0; .25 .25 .25])
zSST = (SSTMETA(14).Data{2}([1 3 5 9],6000:9500)');
zMatSST = zSST;
plot(xMat,zMatSST,'linewidth',3)
ylim([-80 25])
box off
axis off

%We are no longer plotting this.  

% subplot(SP8a)
% plot(1:length(ISISST),ISISST,'-o','Color',GoodCyan,'LineWidth',3,'MarkerFaceColor',GoodCyan)
% box off
% set(gca,'linewidth',4,'fontsize',16,'fontweight','bold')
% title('Inter-Spike Interval SST','FontSize',18)
% SP8a.Title.Position = [8.5 73.7883 0];
% xlabel('ISI Number','FontSize',18)
% % ylabel('ISI(ms)','FontSize',12)
% xlim([0.5 length(ISISST)+0.25])
% ylim([0 70])
% 
% subplot(SP8b)
% boxplot(AdaptIndex,'PlotStyle','compact','Colors',GoodCyan,'Labels',{''});
% h = gca;
% set(h,'linewidth',4,'fontsize',16,'fontweight','bold','yaxislocation','right','yticklabel',[])
% h.Children.Children(6).LineWidth = 2;
% h.Children.Children(5).LineWidth = 8;
% h.XAxis.Visible = 'off';
% ylabel('Adaptation Idex','FontSize',18)
% ylim([-0.025 0.16])
% box off
% 
% subplot(SP9)
% hold on
% errorbar(XSST,AvrSpikesSST,SEMSST,'color','b','LineWidth',4,'Marker','o','MarkerSize',4)
% % scatter([0.42 0.47 0.52],[2.97 4.35 6.1],55,[0 0 1],'d','MarkerFaceColor','b')
% hold off
% box off
% set(gca,'fontsize',10,'fontweight','bold')
% title('Input/Output SST','FontSize',15)
% xlabel('Input Intensity (nA)','FontSize',12)
% ylabel('Spk Count','FontSize',12)
% xlim([XSST(1)*0.95 XSST(end)*1.05])
% ylim([0 18])
% box off

% print '-dtiffn' Figure3
