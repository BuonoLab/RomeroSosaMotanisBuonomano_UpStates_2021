function [Act stimon stimoff dt] = EvokedUp_Fun(GRAPHICS,varargin)
%Outputs:
%Act(time,3) = PV, Ex, SST;
%stimon and stimoff are most for visualization of Opto Stim
%dt is currently fixed in this function.
% CHETA 5 ms on / 10 Off (total of 15 ms)

F = @(x,gain,Thr) gain*max(0,x-Thr); %Threshold Linear activation function for Units ("Input-Output" function);

dt = 0.1;          %theoretical dt, time stepin milliseconds


parser = inputParser;
addParameter(parser,'totSim',3000) %Total Time in Seconds

addParameter(parser,'JEE',5)
addParameter(parser,'JEP',1)
addParameter(parser,'JES',0)
addParameter(parser,'JPE',10)
addParameter(parser,'JPP',0.5)
addParameter(parser,'JPS',0)
addParameter(parser,'JSE',10)
addParameter(parser,'JSP',0)
addParameter(parser,'JSS',0)

addParameter(parser,'thetaE',5) %4.6 [6]   ,  5
addParameter(parser,'thetaP',30)  %25 [36] ,  5(0.36/0.06) = 5*6 = 30
addParameter(parser,'thetaS',15)  %10 [18],   5(0.18/0.06) = 5*3 = 15
addParameter(parser,'gainE',1) %1 [1.24], 1
addParameter(parser,'gainP',2.7) %4 [3.34], 2.7 
addParameter(parser,'gainS',1.6) %1 [1.98]  1.6 (normalized to Py, and gain Py=1 Rocha)

addParameter(parser,'Beta',0) %1 [1.98]
addParameter(parser,'Etau',10) %24 [1.98]
addParameter(parser,'Ptau',4) %8 [1.98]
addParameter(parser,'Stau',6) %15 [1.98]

addParameter(parser,'OUsigmaE',0) %1 [1.98]
addParameter(parser,'OUsigmaP',0) %1 [1.98]
addParameter(parser,'OUsigmaS',0) %1 [1.98]


%% Optogenetics
addParameter(parser,'EvokedOn',500) %1 [1.98]
addParameter(parser,'EvokedDur',25) %1 [1.98]
addParameter(parser,'EvokedAmp',7) %1 [1.98]
addParameter(parser,'OptoAmp_P',0) %1 [1.98]
addParameter(parser,'OptoAmp_S',0) %1 [1.98]
addParameter(parser,'OptoDur',2000);
addParameter(parser,'OptoDelay',250);
addParameter(parser,'seed',2);

OptoTrigThresh = 0.8;
OptoRefact     = 1000/dt; %has to stay above thresh for this long

parse(parser,varargin{:})
p = parser.Results;
totSim = p.totSim/dt;
p.EvokedOn  = p.EvokedOn/dt;
p.EvokedDur = p.EvokedDur/dt;
p.OptoDelay = p.OptoDelay/dt;
p.OptoDur   = p.OptoDur/dt;
p.Etau  = p.Etau/dt;
p.Ptau  = p.Ptau/dt;
p.Stau  = p.Stau/dt;

tauA = 500/dt;    %Time constant of adaptation

% Ornstein–Uhlenbeck Noise
OUtau = 1*dt;
OUmu = 0;

randn('seed',p.seed)
%history
evoked = zeros(1,totSim);
evoked(p.EvokedOn:p.EvokedOn+p.EvokedDur)=p.EvokedAmp;
Act = zeros(totSim,4);
hOpto = zeros(totSim,2);

E = 0; P = 0; S = 0; a = 0;
OUE = 0; OUP = 0; OUS = 0;

OptoStart = -9e9; OptoOff = -9e9;
UpFlag = 0; UpStart = -9e9;
optostim_p = 0; optostim_s = 0;
count = 0;

for t=1:totSim
   count = count+1;
   
   if t>OptoStart && t<=OptoOff
      optostim_p = p.OptoAmp_P;
      optostim_s = p.OptoAmp_S;
   else
      optostim_p = 0;
      optostim_s = 0;
   end
   
   
   OUE = OUE + OUtau*(OUmu-OUE) + p.OUsigmaE*randn; %Ornstein-Uhlenbeck Noise for excitatory unit
   OUP = OUP + OUtau*(OUmu-OUP) + p.OUsigmaP*randn; %Ornstein-Uhlenbeck Noise for inhibitory unit
   OUS = OUS + OUtau*(OUmu-OUS) + p.OUsigmaS*randn; %Ornstein-Uhlenbeck Noise for inhibitory unit
   
   E = E + (-E + F(p.JEE*E - p.JEP*P - p.JES*S- a + OUE + evoked(t),p.gainE,p.thetaE) )/p.Etau;  %Main Equation for Ex firing rate
   P = P + (-P + F(p.JPE*E - p.JPP*P - p.JPS*S + optostim_p + OUP,p.gainP,p.thetaP) )/p.Ptau;
   S = S + (-S + F(p.JSE*E - p.JSP*P - p.JSS*S + optostim_s + OUS,p.gainS,p.thetaS) )/p.Stau;
   
   a = a + (-a + p.Beta*E)/tauA; %ADAPTATION
   
   %Detect Potential Upstate Onset
   if E>OptoTrigThresh && t>OptoOff+OptoRefact
      UpFlag = UpFlag+1;
   else
      UpFlag = 0;
   end
   %Set OptoStim Duraton of Up lasted OptoTrigDur
   if UpFlag == p.OptoDelay
      OptoStart = t;
      OptoOff   = t+p.OptoDur;
   end
   
   Act(t,:) = [P E S a];
   hOpto(t,:) = [optostim_p optostim_s];  %optical stimulation history
   
end

stimon = find(diff(abs(sum(hOpto,2)))>0.001);
stimoff = find(diff(abs(sum(hOpto,2)))<-0.001);
if length(stimoff) < length(stimon)
   stimoff(end+1) = totSim/dt;
end

if GRAPHICS
   [UpFreq UpDur] = RasterFindUpStates(Act(:,2)',1,dt);
   MeanUpDur = mean(UpDur);
   
   
   cla
   colormap = [1 0 0; 0 0.5 0; 0 .75 .75; 0.5 0.5 0.5];
   set(gca,'colororder',colormap);
   hold on
   
   %[UpFreq UpDur] = RasterFindUpStates(hR(:,2)',1);
   taxis = (dt:dt:totSim*dt)/1000;
   plot(taxis,Act(:,1:3),'linewidth',2,'linestyle','-')
   %   plot(taxis,Act(:,4),'linewidth',2,'color',[0.5 0.5 0.5])
   str = sprintf('Red = P; Green=Ex; Blue=S');
   %str = sprintf('Mean UpDur=%6.4f NumUp=%4d',mean(UpDur)*dt,length(UpDur));
   title(str)
   xlabel('Time (s)','fontweight','bold','fontsize',14);
   ylabel('Hz (Red,Green)','fontweight','bold','fontsize',14)
   
   if length(stimon)>length(stimoff)
      stimoff(length(stimon)) = totSim;
   end
   if sum(hOpto(:))>0
      lasercolor = [0 0 1];
   elseif sum(hOpto(:))<0
      lasercolor = 'y';
   end
   
   ylim([0 40]);
   xlim([0 taxis(end)]);
   xl = xlim;
   yl = ylim;
   
   % PLOT OPTOBARS
   
   for i=1:length(stimon)
      hf = fill([stimon(i) stimoff(i) stimoff(i) stimon(i) stimon(i)]*dt/1000,[yl(1) yl(1) yl(2) yl(2) yl(1)],lasercolor);
      set(hf,'edgecolor','none','facealpha',0.25)
   end
   
   xlim(xl)
   drawnow
end

