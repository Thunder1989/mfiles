% ZC number of zero crossings in x
% [n] = zc(x) calculates the mean time scale using generalized zero
% crossing GZC
% type is either 'mf', 'wmf' or 'mp' for mean frequency, weighted mean
% frequency or mean period

function [res] = gzc(data,aggBin,type)
 %%
 
weight = 1;
if strcmp(type,'wmf')
    weight=0;
end

data = interp1(data,[1:1/aggBin:length(data)]);
s=sign(data);
t=filter([1 1],1,s);
t(t==1) = 0;

zcIndex = find(t==0);
timeScales = zeros(7,length(zcIndex)*2-9);
i=1;

if ~isempty(timeScales)
    % for each zero crossing points x
    for ind = 3:1:length(zcIndex)-2
        x = zcIndex(ind);

        % local min and local max preceding x
        [C,localMax] = max(data(zcIndex(ind-2):zcIndex(ind)));
        [C,localMin] = min(data(zcIndex(ind-2):zcIndex(ind)));

        xm3 = min([localMax localMin])+zcIndex(ind-2)-1;
        xm2 = zcIndex(ind-1);
        xm1 = max([localMax localMin])+zcIndex(ind-2)-1;


        % local min and max following x
        [C,localMax] = max(data(zcIndex(ind):zcIndex(ind+2)));
        [C,localMin] = min(data(zcIndex(ind):zcIndex(ind+2)));

        xp1 = min([localMax localMin])+zcIndex(ind)-1;
        xp2 = zcIndex(ind+1);
        xp3 = max([localMax localMin])+zcIndex(ind)-1;
        xp4 = zcIndex(ind+2);

        timeScales(1,i) = (1+weight*3)*(xp1-x);
        timeScales(2,i) = (1+weight)*(xp1-xm1);
        timeScales(3,i) = (1+weight)*(xp2-x);
        timeScales(4,i) = xp1-xm3;
        timeScales(5,i) = xp2-xm2;
        timeScales(6,i) = xp3-xm1;
        timeScales(7,i) = xp4-x;

        if ind~=length(zcIndex)-2
            i=1+i;

            xm3 = xm2;
            xm2 = xm1;
            xm1 = x;
            x   = xp1;
            xp1 = xp2;
            xp2 = xp3;
            xp3 = xp4;

             % local min or max following x+2
            [C,localMax] = max(data(zcIndex(ind+1):zcIndex(ind+3)));
            [C,localMin] = min(data(zcIndex(ind+1):zcIndex(ind+3)));

            xp4 = max([localMax localMin])+zcIndex(ind+1)-1;       

            timeScales(1,i) = (1+weight*3)*(xp1-x);
            timeScales(2,i) = (1+weight)*(xp1-xm1);
            timeScales(3,i) = (1+weight)*(xp2-x);
            timeScales(4,i) = xp1-xm3;
            timeScales(5,i) = xp2-xm2;
            timeScales(6,i) = xp3-xm1;
            timeScales(7,i) = xp4-x;
        end

        i=1+i;
    end

    if strcmp(type,'wmf')
        % return the weighted mean frequency
        res = 1/12*1/mean(sum(timeScales));
    elseif strcmp(type,'mf')
        % return the mean frequency
        res = 1/mean(mean(timeScales));
    else
        %return the mean period
        res = mean(mean(timeScales));
    end
else
    res = Inf; %The period is longer than the data: return infinity
    
end
% figure(110)
% plot(data);
% 
% figure(120)
% plot(median(timeScales));
% hold on
% plot([1 length(timeScales)],[mean(median(timeScales)) mean(median(timeScales))],'g-');
% plot(mean(timeScales),'k');    % non-weigthed mean (see Huang's patent on generalized zero crossing)
% plot([1 length(timeScales)],[mean(mean(timeScales)) mean(mean(timeScales))],'r-');
% pause