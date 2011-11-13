function [a, b, meanDiff] = calcPara(datafile, DEFile)
% This function calculates the mean and variance for UDIFF as it will be
% required to calculate synthetic data.

load(datafile);
load(DEFile);

% First calculate the data control and datanoncontrol.

[row, col] = size(data);
totalSample = col;

timZeroCol = 3:3:totalSample;
timeOneCol = 2:3:totalSample;
timeTwoCol = 1:3:totalSample;

dataZero = data(:,timZeroCol);
dataOne = data(:,timeOneCol);
dataTwo = data(:, timeTwoCol);


% This data is taken before the radiation is applied.
dataControl = dataZero;

% For a gene select data for dataNonControl from the time point that gives a better distance
% from the dataControl.

%dataNonControl = dataOne;


for i=1:row
    % For each of the gene
    avgCtrl = mean(dataZero(i,:));
    
    can1 = mean(dataOne(i,:));
    can2 = mean(dataTwo(i,:));
    
    diff1 = abs(can1 - avgCtrl);
    diff2 = abs(can2 - avgCtrl);
    
    if diff1 > diff2
        dataNonControl(i,:) = dataOne(i,:);
    else
        dataNonControl(i,:) = dataTwo(i,:);
    end
end


for i=1:length(DE)
    gene = DE(i);
    meanDC = mean(dataControl(gene,:));
    meanNDC = mean(dataNonControl(gene, :));
    
    meanDiff(i) = abs(meanDC - meanNDC);
    
    % Pick 
end


% Now calculate the mean and standard deviation of meanDiff

[a, b] = gamfit(meanDiff);

%meanVal = mean(meanDiff);
%sigmaVal = std(meanDiff);

end