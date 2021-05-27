function [UpFreq UpDur] = RasterFindUpStates(Data,threshold,dt);
%data is a 2D raster (cell x time)
%threshold is in Z units

if nargin==1
   threshold = 0.1;
end

if nargin<3
   dt = 1;
end

medfiltwindow = 100/dt; %ms
maxIEI = 50/dt; %max Inter-Event Interval in ms (below this Interval Events will be the same UpState

[numCells maxt] = size(Data);
maxtseconds  = (maxt/1000)*dt;

%% Find Event onsets
%ztraces = zscore(Data')';

%base = median(ztraces')';
%ztraces = ztraces - base;

traces = mean(Data,1);

%meanz    = mean(ztraces,1);      % mean of the zscore of all traces
smoothz  = medfilt1(traces,medfiltwindow,'truncate');    % smooth mean of zscores
smoothz  = smoothz - min(smoothz);

above        = find(smoothz>threshold);
eventindexes = find(diff([-9e9 above])>1);
eventindexes = find(diff([-9e9 above])>maxIEI);

if ~isempty(above)
   off      = above([eventindexes(2:length(eventindexes))-1 length(above)]);
   Onsets   = above(eventindexes);
else
   off = [];
   Onsets = [];
end

if length(off)>0 && off(end) == maxt
   off = off(1:end-1);
   Onsets = Onsets(1:end-1);
end

UpNumber = length(Onsets);
UpFreq = UpNumber/maxtseconds;

UpDur  = off-Onsets;
UpDur  = UpDur*dt;
Onsets = Onsets*dt;
off    = off*dt;

if (0)
   figure
   plot((1:length(smoothz))*dt,smoothz)
   line([get(gca,'Xlim')],[threshold threshold],'linestyle',':','color','k','linewidth',1);
   ylim = get(gca,'Ylim'); set(gca,'Ylim',[-1 ylim(2)]); ylim = get(gca,'Ylim');
   line([Onsets; Onsets],repmat([ylim(1) ylim(2)]',1,size(Onsets,1)),'linestyle',':','color','g','linewidth',1);
   line([off; off],repmat([ylim(1) ylim(2)]',1,size(Onsets,1)),'linestyle',':','color','r','linewidth',1);
end