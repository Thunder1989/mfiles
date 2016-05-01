file = dir('/Users/hdz_1989/Documents/MATLAB/GPS_Data/*.xlsx');
figure
for n = 1:size(file,1)
    str = ['/Users/hdz_1989/Documents/MATLAB/GPS_Data/', file(n).name]
    %[data, txt] = xlsread(str,'Home');
    data_cleaning_loose;
    
    subplot(4,2,n)
    h = bar(energy_penalty(1:step+2,:),'stack');
    set(h,'BarWidth',0.6); % The bars will now touch each other
    set(gca,'YGrid','on')
    %set(gca,'GridLineStyle','+')
    %set(gca,'XTicklabel','LW>=sD&&LH>=sD|LW>=sD&&LH<=sD|LW<=sD&&LH>=sD|LW<=sD&&LH<=sD')
    set(gca,'XTicklabel','5|10|15|20|25|30|baseline|preheat')
    set(get(gca,'YLabel'),'String','Energy Consumption/kWh')
    lh = legend('Maintain','React');
    set(lh,'Location','SouthEast','Orientation','vertical')
    colormap(summer) 
end