% RUN Figure 3 first


%add zeros column
StepDur = .25; %sec

%% Init Graphics
figure
set(gcf,'position',[80  60  1160  380],'color','w')
SP(1) = subplot('Position',[0.052 0.175 0.24 0.72]);
SP(2) = subplot('Position',[0.37 0.175 0.24 0.72]);
SP(3) = subplot('Position',[0.705 0.175 0.24 0.72]);


%% Main Plotting
for celltype = 1:3
   
   switch celltype
      case 1
         Spikes = [zeros(size(SpikesPYR,1),1) SpikesPYR];
         [numCells maxTraces] = size(SpikesPYR);
         color = [0 0.7 0];
         text = 'Pyr';
         ymax = 50;
      case 2
         Spikes = [zeros(size(SpikesPV,1),1) SpikesPV];
         [numCells maxTraces] = size(SpikesPV);
         color = 'r';
         text = 'PV';
         ymax = 100;
      case 3
         Spikes = [zeros(size(SpikesSST,1),1) SpikesSST];
         [numCells maxTraces] = size(SpikesSST);
         color = [0 0.9 0.9];
         text = 'SST';
         ymax = 80;
   end
   INTEN = [0 0.05:0.05:0.05*(maxTraces)];
   Inten{celltype} = INTEN;
   SPIKES{celltype} = Spikes;
   Spks = nan(numCells,length(INTEN));
   
   count = 1; Accept = [];
   for i=1:numCells
      num = sum(Spikes(i,:));
      if num>1
         Accept = [Accept i];
         Spks(count,:) = Spikes(i,:);
         count = count+1;
      end
   end
   
   Spks    = Spks(1:count-1,:);
   Spks    = Spks/StepDur; %Hz
   SpkMean = nanmean(Spks,1);
   SpkSEM  =  nanstd(Spks)/sqrt(size(Spks,1));
   
   subplot(SP(celltype))
   hold on
   labelfontsize = 15;
   xlabel('Intensity (nA)','fontweight','bold','fontsize',labelfontsize)
   if celltype == 1
      ylabel('Firing rate','fontweight','bold','fontsize',labelfontsize)
   end
   box off
    if celltype ~= 3
        plot(INTEN,Spks','color',[0.8 0.8 0.8])
    elseif celltype ==3
        plot(INTEN,Spks([1:4 6:end],:)','color',[0.8 0.8 0.8])
    end
   hold on
   
   errorbar(INTEN,SpkMean,SpkSEM,'linewidth',5,'color','k')
   %set(gca,'ylim',[0 50])
   
   %% FIT
   F = @(param,x) param(1)*max(0,x-param(2)); %Threshold Linear activation function for Units ("Input-Output" function);
   X = INTEN;
   Y = SpkMean;
   
   param0 = [0, 0 ];
   param  = lsqcurvefit(F,param0,X,Y);
   Param(celltype,:) = param;
   
   %% Plot Fit of Mean Excitability
   X = X(1):0.01:X(end);
   Yhat = F(param,X);
   plot(X,Yhat,'linewidth',3,'color',color);
   
   set(gca,'linewidth',3,'fontweight','bold','fontsize',16);
   set(gca,'xlim',[0 INTEN(end)],'ylim',[0 ymax]);
   str = sprintf('\\fontsize{18}%s: \\fontsize{16}g=%3d, \\theta=%3.2f',text,round(Param(celltype,1)),Param(celltype,2));
   title(str)
   
   set(gca,'ylim',[0 75])
   set(gca,'xlim',[0 0.55])
   
end
% 
% set(gcf,'paperpositionmode','auto','color','w')
% print -djpeg100 ExcitFits.jpg -r300

%% Stats 

ModelSPIKES     = [reshape(SPIKES{1}',[],1); reshape(SPIKES{2}',[],1); reshape(SPIKES{3}',[],1)];
ModelINTEN      = [reshape(repmat(Inten{1},size(SPIKES{1},1),1)',[],1); reshape(repmat(Inten{2},size(SPIKES{2},1),1)',[],1); reshape(repmat(Inten{3},size(SPIKES{3},1),1)',[],1)];
ModelCELLTYPE   = [reshape(repmat({'PYR'},size(SPIKES{1},1),size(SPIKES{1},2))',[],1); reshape(repmat({'PV'},size(SPIKES{2},1),size(SPIKES{2},2))',[],1); reshape(repmat({'SST'},size(SPIKES{3},1),size(SPIKES{3},2))',[],1)];
ModelUNIT       = [reshape(repmat((1:size(SPIKES{1},1))',1,size(SPIKES{1},2))',[],1); reshape(repmat((size(SPIKES{1},1)+1:size(SPIKES{1},1)+size(SPIKES{2},1))',1,size(SPIKES{2},2))',[],1); reshape(repmat((size(SPIKES{1},1)+size(SPIKES{2},1)+1:size(SPIKES{1},1)+size(SPIKES{2},1)+size(SPIKES{3},1))',1,size(SPIKES{3},2))',[],1)];

DataTable       = table(ModelSPIKES,ModelINTEN,ModelCELLTYPE,ModelUNIT); 
ModelFormula    = 'ModelSPIKES~1+ModelINTEN*ModelCELLTYPE+(1+ModelINTEN|ModelUNIT)';

%With a Linear Mixed Effects Model (We assume that the data is normalized distributed)
LinearMixedModel  = fitlme(DataTable,ModelFormula);
disp('Linear Mixed Effects Model')
anova(LinearMixedModel)

%With a Generalized Linear Mixed Effects Model (Doesn't assume normalization)
GeneralizedLinearMixedModel  = fitglme(DataTable,ModelFormula);
disp('Generalized Linear Mixed Effects Model')
anova(GeneralizedLinearMixedModel)

