%load SearchOut_REF3
% load SearchOut_REF4
%load SearchOut_REF5
%load SearchOut_REF6
% cc
close all
clear

F = @(x,gain,Thr) gain*max(0,x-Thr);

filename = 'EvokedParamSearch_FullConnect_REF1_Lower&UpperBounds';
% filename = 'EvokedParamSearch_FullConnect_REF1_UpperBounds';
% filename = 'JL_ReplicatingResultsV2'
% filename = 'ExciteLowerBoundsReviewer'

if (0)
   
   load(filename)
   %load SearchOut_REF5_13^9_tau10_4_6
   %load SearchOut_REF2_13^9_tau15_5_15
   %load SearchOut_REF3_13^9_tau15_5_5
   
   EvokedOn = 500;
   OptoDelay = 250;
   dt = 0.1;
   
   figure
   set(gcf,'position',[  100    300    1100   340])
   SP1 = subplot(2,2,1);
   SP2 = subplot(2,2,2);
   SP3 = subplot(2,2,3);
   SP4 = subplot(2,2,4);
   
   
   %% GRAB PARAMS WITH +/- 25% and stable
   THR = 0.25;
   if (0)
      stdthr = STD<0.1;             %cutoff for oscilating/ringing networks
      ssP = abs(SS(:,1)-TARGET(1))<TARGET(1)*THR;
      ssE = abs(SS(:,2)-TARGET(2))<TARGET(2)*THR;
      ssS = abs(SS(:,3)-TARGET(3))<TARGET(3)*THR;
      ssthr  = ssP.*ssE.*ssS;
      accept = stdthr.*ssthr';
      AcceptInd = find(accept);
   end
   
   pcount = 0; pcountS = 0; PotWP = []; PotWS = [];
   MaxEig = zeros(1,length(AcceptInd));
   META_SS       = zeros(length(AcceptInd),3);
   META_OptoAct1 = zeros(length(AcceptInd),3);
   META_OptoAct2 = zeros(length(AcceptInd),3);
   GRAPHICS = 1;
   
   for i = 1:length(AcceptInd) %1:250 %2047:2049
   %parfor i = 1:length(AcceptInd) %1:250
      
      best = AcceptInd(i);
      
      jee=WEIGHTS(best,1);
      jep=WEIGHTS(best,2);
      jes=WEIGHTS(best,3);
      jpe=WEIGHTS(best,4);
      jpp=WEIGHTS(best,5);
      jps=WEIGHTS(best,6);
      jse=WEIGHTS(best,7);
      jsp=WEIGHTS(best,8);
      jss=WEIGHTS(best,9);
      
      fprintf('%3d: JEE = %4.2f',best,jee);
      fprintf('   JEP = %4.2f',jep);
      fprintf('   JES = %4.2f',jes);
      fprintf('   JPE = %4.2f',jpe);
      fprintf('   JPP = %4.2f',jpp);
      fprintf('   JPS = %4.2f',jps);
      fprintf('   JSE = %4.2f',jse);
      fprintf('   JSP = %4.2f',jsp);
      fprintf('   JSS = %4.2f\n',jss);
      
      if GRAPHICS, subplot(SP1), end
      
      %      [Act1, ~, ~, dt]= EvokedUp_Fun(1,...
      %                            'JEE',jee,...
      %                            'JEP',jep,...
      %                            'JES',jes,...
      %                            'JPE',jpe,...
      %                            'JPP',jpp,...
      %                            'JPS',jps,...
      %                            'JSE',jse,...
      %                            'JSS',jss,...
      %                            'totSim',2000);
      
      Act1 = EvokedUp_FullConnect_Fun(GRAPHICS,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSP',jsp,'JSS',jss,...
         'EvokedOn',EvokedOn,'OptoDelay',OptoDelay,'OptoAmp_P',5,'OptoDur',250,'totSim',1500);
      
      
      if GRAPHICS
         subplot(SP3),
         ylabel('PV','fontweight','bold')
      end
      
      Act2 = EvokedUp_FullConnect_Fun(GRAPHICS,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSP',jsp,'JSS',jss,...
         'EvokedOn',EvokedOn,'OptoDelay',OptoDelay,'OptoAmp_S',5,'OptoDur',250,'totSim',1500);
      ylabel('SST','fontweight','bold')
      
      %    subplot(SP2)
      %    Act = EvokedUp_Fun(1,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSS',jss,...
      %       'OptoAmp_P',-5);
      %
      %    subplot(SP4)
      %    Act = EvokedUp_Fun(1,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSS',jss,...
      %       'OptoAmp_S',-5);
      
      ss = mean(Act1(end-500:end,1:3)); %SteadyState
      mse = sum((TARGET(1:3)-ss(1:3)).^2);
      stdev = std( Act1((EvokedOn+OptoDelay-100)/dt:(EvokedOn+OptoDelay-1)/dt,2) );
      optoact1 = mean( Act1( (EvokedOn+OptoDelay+100:EvokedOn+OptoDelay+200)/dt,1:3) );
      optoact2 = mean( Act2( (EvokedOn+OptoDelay+100:EvokedOn+OptoDelay+200)/dt,1:3) );
      
      if mean(Act1(end-500:end,2))<0.1
         PotWP = [ PotWP best];
      end
      if mean(Act2(end-500:end,2))<0.1
         PotWS = [ PotWS best];
      end
      
      % Analytical Steady State
      
      %JEE_p = jee - 1/gainE;
      %JII_p = jpp + 1/gainP;
      %M = jep*jpe - (JEE_p-Beta)*(JII_p);
      
      %SSE =  1/M*(jep*thetaP - JII_p*thetaE);
      %SSI =  1/M*((JEE_p-Beta)*thetaP - jpe*thetaE);
      
      E = ss(2); P = ss(1); S = ss(3);
      gE = 1; gP = 2.7; gS = 1.6;
      sse = (-E + F(E*jee - P*jep - S*jes,gE,5 ) )/10;
      ssp = (-P + F(E*jpe- P*jpp - S*jps,gP,30 ))/4;
      sss = (-S + F(E*jse- P*0 - S*jss,gS,15 ))/6;
      jeePrime = jee-1/gE;
      jppPrime = jpp+1/gP;
      jssPrime = jss+1/gS;
      
      W = [jee -jep -jes; jpe -jpp -jps; jse  -jsp -jss];
      THETA = diag([5 30 15]);
      GAIN = diag([1 2.7 1.6]);
      T = diag([10 4 6]);
      [V, D] = eig(T^-1*(-eye(size(W,1))+GAIN*W));
      MaxEig(i) = max(diag(real(D)));
      META_OptoAct1(i,:) = optoact1;
      META_OptoAct2(i,:) = optoact2;
      META_SS(i,:)       = ss;
      META_PredictPrdx(i,1) = jeePrime*jssPrime - jes*jse;
      META_PredictPrdx(i,2) = jeePrime*jppPrime - jep*jpe;
      META_ActualPrdx(i,1) = ss(1)-optoact1(1);
      META_ActualPrdx(i,2) = ss(3)-optoact2(3);
      
      fprintf('%6d/%6d| SS = %4.2f %4.2f %4.2f | STD=%4.2f | Max Eig = %6.3f\n',i,length(AcceptInd),ss,stdev,MaxEig(i));   %fprintf('SSE/I = %4.2f/%4.2f/%4.2f | Analytical SS = %4.2f/%4.2f\n',SS(2),SS(1),SS(3),SSE,SSI);
      fprintf('%6d| delta_SS(empricial-analytical)= %4.2f %4.2f %4.2f\n',i,sse,ssp,sss);   %fprintf('SSE/I = %4.2f/%4.2f/%4.2f | Analytical SS = %4.2f/%4.2f\n',SS(2),SS(1),SS(3),SSE,SSI);
      fprintf('Predicted/Actual P Paradoxical Effect = %6.3f/%6.3f\n',META_PredictPrdx(i,1),META_ActualPrdx(i,1))
      
      if GRAPHICS
         waitforbuttonpress
         drawnow
      end
      
      fprintf('\n')

      
   end
      %% PERCENT OF WEIGHT SET that turned off Up-states
      PerPVW = size(PotWP,2)/length(AcceptInd)
      PerSSTW = size(PotWS,2)/length(AcceptInd)
   
end
%% MINE THROUGH THE PotW Weights to see if any exhibit Longer UpStates with

if 0
   
   figure
   set(gcf,'position',[  200    400    1100   340])
   SP1 = subplot(2,2,1);
   SP2 = subplot(2,2,2);
   SP3 = subplot(2,2,3);
   SP4 = subplot(2,2,4);
   
   clear UpDur1 UpDur3 UpDur4
   for best = 1:size(PotWP,1)
      
      jee=PotWP(best,1);
      jep=PotWP(best,2);
      jes=PotWP(best,3);
      jpe=PotWP(best,4);
      jpp=PotWP(best,5);
      jps=PotWP(best,6);
      jse=PotWP(best,7);
      jsp=PotWP(best,8);
      jss=PotWP(best,9);
      jee_p = jee
      
      fprintf('   JEE = %4.2f',jee);
      fprintf('   JEP = %4.2f',jep);
      fprintf('   JES = %4.2f',jes);
      fprintf('   JPE = %4.2f',jpe);
      fprintf('   JPP = %4.2f',jpp);
      fprintf('   JPS = %4.2f',jps);
      fprintf('   JSE = %4.2f',jse);
      fprintf('   JSP = %4.2f',jsp);
      fprintf('   JSS = %4.2f\n',jss);
      
      
      OUsigmaE = 1.4; totSim = 5000;
      OUsigmaP = 4;
      subplot(SP1)
      Act = EvokedUp_FullConnect_Fun(1,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,'JSP',jsp,'JSS',jss,...
         'OptoAmp_P',0,'OUsigmaE',OUsigmaE,'OUsigmaP',OUsigmaP,'totSim',totSim,'EvokedAmp',0);
      [UpFreq updur] = RasterFindUpStates(Act(:,2)',1,dt);
      UpDur1(best) = mean(updur);
      
      
      
      subplot(SP2)
      Act3 = EvokedUp_FullConnect_Fun(1,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,jse,'JSP','JSS',jss,...
         'OptoAmp_P',-10,'OUsigmaE',OUsigmaE,'OUsigmaP',OUsigmaP,'totSim',totSim,'EvokedAmp',0);
      [UpFreq updur] = RasterFindUpStates(Act3(:,2)',1,dt);
      UpDur3(best) = mean(updur);
      
      subplot(SP4)
      Act4 = EvokedUp_FullConnect_Fun(1,'JEE',jee,'JEP',jep,'JES',jes,'JPE',jpe,'JPP',jpp,'JPS',jps,'JSE',jse,jse,'JSP','JSS',jss,...
         'OptoAmp_S',-10,'OUsigmaE',OUsigmaE,'OUsigmaP',OUsigmaP,'totSim',totSim,'EvokedAmp',0);
      [UpFreq updur] = RasterFindUpStates(Act4(:,2)',1,dt);
      UpDur4(best) = mean(updur);
      
      SS = mean(Act(end-1000:end,1:3)); %SteadyState
      MSE = sum((TARGET(1:3)-SS(1:3)).^2);
      
      
      
      % Analytical Steady State
      
      
      fprintf('%2d | Dur: %4.2f %4.2f %4.2f\n',best,UpDur1(best),UpDur3(best),UpDur4(best));
      fprintf('\n')
      
      waitforbuttonpress
      drawnow
      
      
   end
   
end






%% WEIGHT DISTRIBUTIONS
if (1)
   
   load(filename)
   %load SearchOut_REF5_13^9_tau10_4_6
   %load SearchOut_REF2_13^9_tau15_5_15
   %load SearchOut_REF3_13^9_tau15_5_5
   
   
   %% GRAB PARAMS WITH +/- 20% and stable
   stdthr = STD<0.1;
   ssP = abs(SS(:,1)-TARGET(1))<TARGET(1)*0.25;
   ssE = abs(SS(:,2)-TARGET(2))<TARGET(2)*0.25;
   ssS = abs(SS(:,3)-TARGET(3))<TARGET(3)*0.25;
   ssthr  = ssP.*ssE.*ssS;
   accept = stdthr.*ssthr';
   AcceptInd = find(accept);
   
   
   fontsize = 12;
   W = WEIGHTS(AcceptInd,:);
   %   W = PotWP;
   %   W = PotWS;
   
   
   IP = W(:,2).*W(:,4);
   IS = W(:,3).*W(:,7);
   
   
   %% Net P Inhibition - Net S Inhibition
   
   figure, set(gcf,'color','w','position',[200  400  460  420])
   binsize = 2; %This doesn't really matter since the actual binsizes are hardcoded on the histogramline...
   h1 = histogram(IP-IS,[min(IP-IS)-binsize/2:4:max(IP-IS)+binsize/2],'facecolor','b');
   
   
   set(gca,'color','w','linewidth',2,'fontweight','bold','ycolor','k','fontsize',fontsize-2,'xcolor','k')
   set(h1,'facecolor',[0.5 0.5 0.5]);
   ylabel('# Weight Sets','fontweight','bold','fontsize',fontsize+6)
   xlabel(['{\color{red}W_E_P*W_P_E} - {\color[rgb]{0 .9 .9}W_E_S*W_S_E}'],'fontweight','bold','fontsize',fontsize+4)
   line([0 0 ],ylim,'color',[0.2 0.2 0.2],'linestyle',':','linewidth',2)
   box off
   
   %set(gcf,'inverthardcopy','off')
   set(gcf,'paperpositionmode','auto','color','w')
   print('Dist1','-djpeg','-r300')
   
   %% NET PINHIBITION - NET P EXCITATION
   fontsize = 12;
   
   IP = W(:,2).*W(:,4);
   EP = W(:,1).*W(:,5);
   
   figure, set(gcf,'color','w','position',[220  400  460  420])
   h1 = histogram(IP-EP,[min(IP-EP)-binsize/2:4:max(IP-EP)+binsize/2],'facecolor','b');
   
   
   set(gca,'color','w','linewidth',2,'fontweight','bold','ycolor','k','fontsize',fontsize-2,'xcolor','k')
   set(h1,'facecolor',[0.5 0.5 0.5]);
   ylabel('# Weight Sets','fontweight','bold','fontsize',fontsize+6)
   xlabel('{\color{red}W_E_P*W_P_E} - {\color[rgb]{0 .7 0}W_E_E*W_P_P}','fontweight','bold','fontsize',fontsize+4)
   line([0 0 ],ylim,'color',[0.2 0.2 0.2],'linestyle',':','linewidth',2)
   box off
   
   %set(gcf,'inverthardcopy','off')
   set(gcf,'paperpositionmode','auto','color','w')
   print('Dist2','-djpeg','-r300')
   
   %% NET S INHIBITION - NET S EXCITATION
   
   IS = W(:,3).*W(:,7);
   ES = W(:,1).*W(:,9);
   
   figure, set(gcf,'color','w','position',[240  400  460  420])
   h1 = histogram(IS-ES,[min(IS-ES)-binsize/2:4:max(IS-ES)+binsize/2],'facecolor','b');
   
   set(gca,'color','w','linewidth',2,'fontweight','bold','ycolor','k','fontsize',fontsize-2,'xcolor','k')
   set(h1,'facecolor',[0.5 0.5 0.5])
   ylabel('# Weight Sets','fontweight','bold','fontsize',fontsize+6)
   xlabel('{\color[rgb]{0 0.9 0.9}W_E_S*W_S_E} - {\color[rgb]{0 0.7 0}W_E_E*W_S_S}','fontweight','bold','fontsize',fontsize+4)
   line([0 0 ],ylim,'color',[0.2 0.2 0.2],'linestyle',':','linewidth',2)
   box off
   
   %set(gcf,'inverthardcopy','off')
   set(gcf,'paperpositionmode','auto','color','w')
   print('Dist3','-djpeg','-r300')
   
end


%% CLUSTERING
if (0)
   %WC = W./max(W);
   WC = W;
   %WC = [W(:,1) W(:,2).*W(:,4) W(:,1).*W(:,5) W(:,3).*W(:,7) W(:,1).*W(:,9) W(:,6) W(:,8)];
   [kclust] = kmeans(WC,2,'distance','corr','replicates',250);
   [val ind]=sort(kclust);
   imagesc(WC(ind,:))
   
   E = evalclusters(WC,'kmeans','Silhouette','Distance','corr','klist',[1:5]);
   
end


