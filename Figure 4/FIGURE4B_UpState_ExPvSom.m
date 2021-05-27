close all
clear all

%load SearchOut_SearchFullConnect_REF1_46^6 [
load EvokedParamSearch_FullConnect_REF1_Lower&UpperBounds

dt = 0.1;          %theoretical dt, time stepin milliseconds
totSim = 16000;
%totSim = 20000;
OUsigma = 1;
Beta = 1; %0.8
taxis = (dt:dt:totSim)/1000;

close all

r = [1 0.1 0.1]; g = [0.25 1 0.25]; c = [0 0.95 0.95];

%% LOOP THROUGH W's w/ Lowest MSE's

THR = 0.25;
if (1)
   stdthr = STD<0.1;             %cutoff for oscilating/ringing networks
   ssP = abs(SS(:,1)-TARGET(1))<TARGET(1)*THR;
   ssE = abs(SS(:,2)-TARGET(2))<TARGET(2)*THR;
   ssS = abs(SS(:,3)-TARGET(3))<TARGET(3)*THR;
   ssthr  = ssP.*ssE.*ssS;
   accept = stdthr.*ssthr';
   AcceptInd = find(accept);
end

% FIND WEIGHT SET CLOSEST TO THE SET POINT
AccW = WEIGHTS(AcceptInd,:);
CentroidW=mean(AccW);
CentroidDist = pdist2(CentroidW,AccW);
[v Prototype] = min(CentroidDist);
[v sanitycheck]=min(sum(pdist2(AccW,AccW)));


%% PLOT PRETTY FIGURE SPONTANEOUS (WOB)
for best = Prototype
      
   ind = AcceptInd(best);
      
   jee=WEIGHTS(ind,1); jep=WEIGHTS(ind,2); jes=WEIGHTS(ind,3);
   jpe=WEIGHTS(ind,4); jpp=WEIGHTS(ind,5); jps=WEIGHTS(ind,6);
   jse=WEIGHTS(ind,7); jsp=WEIGHTS(ind,8); jss=WEIGHTS(ind,9);

   fprintf('%3d: JEE = %4.2f',best,jee);
   fprintf('   JEP = %4.2f',jep);
   fprintf('   JES = %4.2f',jes);
   fprintf('   JPE = %4.2f',jpe);
   fprintf('   JPP = %4.2f',jpp);
   fprintf('   JPS = %4.2f',jps);
   fprintf('   JSE = %4.2f',jse);
   fprintf('   JSP = %4.2f',jsp);
   fprintf('   JSS = %4.2f\n',jss);
   
   %    Act = EvokedUp_FunNew(0,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSS',jss,...
   %       'OUsigmaE',0.9*OUsigma,'OUsigmaP',3*OUsigma,'OUsigmaS',1*OUsigma,'totSim',totSim, ...
   %       'EvokedAmp',0,'Beta',Beta);
   
   ESigmaFactor = 1.1;
   ISigmaFactor = 5;
   [Act stimon stimoff]= EvokedUp_FullConnect_Fun(0,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSP',jsp,'JSS',jss,...
      'OUsigmaE',ESigmaFactor*OUsigma,'OUsigmaP',ISigmaFactor*OUsigma,'OUsigmaS',ISigmaFactor*OUsigma,'totSim',totSim, ...
      'EvokedAmp',0,'Beta',Beta*0,'seed',57); %20 34 68 63 57
   %The line above is where we randomly went through many seeds to find the one we used for the figure
   
   hR = medfilt1(Act(:,:),50/dt,'truncate');
   
   yshiftP = 20;
   yshiftS = 20;
   
   figure
   set(gcf,'position',[80  270  560.0000  290],'color','w')
   fontsize = 20;
   plot(taxis,hR(:,2),'color',[g],'linewidth',2)
   hold on
   plot(taxis,hR(:,1)-yshiftP,'color',[r 0.9],'linewidth',2)
   plot(taxis,hR(:,3)-yshiftP-yshiftS,'color',[c],'linewidth',2)
   ylabel('Firing rate (Hz)','fontweight','bold','fontsize',fontsize)
   set(gca,'color','w','linewidth',2,'fontweight','bold','ycolor','k','fontsize',14,'xcolor','k')
   xlabel('Time (s)','fontweight','bold','fontsize',fontsize)
   box off
   str = sprintf('UpDown_BOW%d',1);
   
   set(gca,'ylim',ylim-[5 0])
   set(gca,'ytick',[-yshiftP-yshiftS -yshiftP-yshiftS+10 -yshiftP -yshiftP+10  0 10],'yticklabel',{'0', '10', '0', '10','0','10'})
   
   set(gca,'color','w')
   %set(gcf,'inverthardcopy','off')
   set(gcf,'paperpositionmode','auto','color','w')
   print(str,'-djpeg','-r300')
   
end

