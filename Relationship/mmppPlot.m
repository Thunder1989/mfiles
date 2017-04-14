function plotMMPP(L,Z,N,N0,TRUTH,FIG,RANGE)
% plotMMPP(Rate,PEvent,Data,Events,Fig,Range)
% plot data and parameters of the Markov modulated Poisson process
%   Rate   = estimated rate function; 
%   PEvent = estimated probability of event;
%   Data   = observed data; 
%   Events = known event times;
%   Fig    = figure handle / number to use
%   Range  = sub-range of the data to plot (default is full length)

% Copyright (C) 2006 Alexander Ihler; distributable under GPL -- see README.txt

if (~exist('RANGE','var')) RANGE = 1:numel(Z); end;

figure(FIG); 
subplot(3,1,1:2,'replace');
hold on; grid on;
H=plot(RANGE,N(RANGE),'b-');
ylabel('Counts'); %legend('Observed');
% subplot(3,1,2,'replace')
H=plot(RANGE,N0(RANGE),'r-');  
% ylabel('Counts'); legend('Estimated N0');

figure(FIG); subplot(3,1,3,'replace'); hold on;
  S=stem(RANGE,Z(RANGE),'k'); set(S,'MarkerSize',0,'LineWidth',2);
  xlabel('Time'); ylabel('Event');
  if (~isempty(TRUTH)),
    hold on; S=stem(RANGE,-.25*TRUTH(RANGE),'b'); set(S,'MarkerSize',0,'LineWidth',1);
  end;