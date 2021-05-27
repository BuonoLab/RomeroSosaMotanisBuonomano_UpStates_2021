clear all


TARGET = [14 5 17 ]; %Target Hz for PV, E, and SST 5,19,18

%The center of these values come from Jecorg et al 2017
JEE = [2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9];     lee = numel(JEE);
JEP = [0  0.5 1 1.5 2 2.5 3 3.5 4];               lep = numel(JEP); %In the text we say that we go from 0 to 4 but that's not the case in here...
JES = JEP;                                   les = numel(JES);
JPE = [2 4 6 8 10 12 14 16 18];                lpe = numel(JPE); %Why are the increments in this part by 2 where previously they were by 0.5
JPP = [0  1  2  3  4 5];                lpp = numel(JPP); %Same as above in addition, in the text we said we went up to 5 but in reality we are only going up to 4... 
JPS = [0  1  2  3 4 5];                 lps = numel(JPS);
JSE = JPE;                                   lse = numel(JSE);
JSP = JPS;                                   lsp = numel(JSP);
JSS = JPP;                                   lss = numel(JSS);


totalparam = lee*lep*les*lpe*lpp*lps*lse*lsp*lss; %Why are the values in the totalparam here different from the one in the paper (13,505,625 here and 46,294,416)? Figured that this is caused by wrong bounds being set
SS = zeros(totalparam,3);
STD = zeros(1,totalparam);
WEIGHTS = zeros(totalparam,9);

tic
count = 0;
parfor i=1:totalparam
%for i=1:lee*lep*les*lpe*lpp*lps*lse*lss

   [ind9 ind8 ind7 ind6 ind5 ind4 ind3 ind2 ind1] = ind2sub([lss,lsp,lse,lps,lpp,lpe,les,lep,lee],i);
   
   [Act, ~, ~, dt]= EvokedUp_FullConnect_Fun(0,...
                           'JEE',JEE(ind1),...
                           'JEP',JEP(ind2),...
                           'JES',JES(ind3),...
                           'JPE',JPE(ind4),...
                           'JPP',JPP(ind5),...
                           'JPS',JPS(ind6),...
                           'JSE',JSE(ind7),...
                           'JSP',JSP(ind8),...
                           'JSS',JSS(ind9),...
                           'totSim',1500);
                        
                        ss = mean(Act(end-100/dt:end,1:3)); %SteadyState
                        stdev = std( Act(end-100/dt:end,2) );
                        
                        SS(i,:) = ss;
                        STD(i)    = stdev;
                        WEIGHTS(i,:) = [JEE(ind1) JEP(ind2) JES(ind3) JPE(ind4) JPP(ind5) JPS(ind6) JSE(ind7) JSP(ind8) JSS(ind9)];
                        if mod(i,1000)==0
                           fprintf('%6d/%6d:%d:%d:%d:%d:%d:%d:%d:%d %d| %4.2f %4.2f %4.2f \n',i,totalparam,ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,ss(1),ss(2),ss(3));
                        end
end

toc

clear Act*
close all

WEIGHTS = single(WEIGHTS);
SS = single(SS);
STD = single(STD);
save ExciteLowerBoundsReviewer