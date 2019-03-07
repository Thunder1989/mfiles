function plotMMPP(L,Z,N,N0,TRUTH,FIG,RANGE)
% plotMMPP(Rate,PEvent,Data,Events,Fig,Range)
% plot data and parameters of the Markov modulated Poisson process
%   Rate   = estimated rate function; 
%   PEvent = estimated probability of event;
%   Data   = observed data; 
%   Events = known event times;
%   Fig    = figure handle / number to use
%   Range  = sub-range of the data to plot (default is full length)

if (~exist('RANGE','var')), RANGE = 1:numel(Z); end;

figure(FIG); 
grid on; 
subplot(3,1,1:2,'replace');
hold on;
plot(RANGE,N(RANGE),'b-');
% subplot(3,1,2,'replace')
plot(RANGE,N0(RANGE),'r-');  
ylabel('Counts'); %legend('Observed');
% set(gcf, 'YLim', [0 1.5*max(N(RANGE))]);

% figure(FIG); subplot(3,1,3,'replace'); hold on;
% S=stem(RANGE,Z(RANGE),'k'); set(S,'Marker','none','LineWidth',1.5);
% xlabel('Time'); ylabel('Event');
% if (~isempty(TRUTH)),
%     S=stem(RANGE,-.25*TRUTH(RANGE),'b'); set(S,'Marker','none','LineWidth',1);
% end;
