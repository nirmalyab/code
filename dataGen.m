function dataGen(dataFile, maxNeigh, seedArr, jpath, relFile, outFile, intGap, extGap, fixedSigma)
% This function generates data for our experiments. 

% X: contains the parameters for generating data.
% maxNeigh: Contains the maximum number of DE neighbors that a DE gene can
% have. Here neighbors are found using breadth fast search method.
% seedCommon:  The seeds that are common to both the datasets.
% seedA : The seeds that are specific to dataset A.
% seedB: The seeds that are specific to dataset B.
% patientCount: Number of patients.
% intGap and extGap needs to be provided in terms of standard daviation X(2)
% First load the datafile.

load(dataFile);

%[geneCount, junktCount] = size(dataControl);

% Use the same control groups for dataA and dataB.
dataAControl = dataControlAA;
dataBControl = dataControlCA;

dataANonControl = dataNonControlAA;
dataBNonControl = dataNonControlCA;

javaaddpath(jpath);
relProg = RelProg(relFile);
relations = load(relFile);
gRelatedGenes = getRelatedGenes(relations);

len = length(gRelatedGenes);

gGeneStateA = zeros(gRelatedGenes(len),1);
gGeneStateB = zeros(gRelatedGenes(len),1);
gGeneState = zeros(gRelatedGenes(len),1);

% Make all the genes that are related to be 1 (EE).
gGeneState(gRelatedGenes) = 1;
gGeneStateA(gRelatedGenes) = 1;
gGeneStateB(gRelatedGenes) = 1;

% probVal is the initial probability that a gene is DE. Seed is used to
% control the initial genes that that will be made DE or EE. 
inProb = 0.4;

% Make it a cumulative array

countArr = cumsum(seedArr);
[intGeneA, extGeneA, intGeneB, extGeneB, intGeneC, extGeneC] ...
    = initializeGenes(relProg, maxNeigh, sum(seedArr), countArr, ...
    gRelatedGenes, gGeneState, inProb);

% Change the values of gGeneStateA and gGeneStateB based on intGeneA,
% extGeneA, intGeneB, extGeneB, intGeneC, extGeneC. 
GeneADE = [intGeneA, extGeneA, intGeneC, extGeneC];
gGeneStateA(GeneADE) = 2;
GeneBDE = [intGeneB, extGeneB, intGeneC, extGeneC];
gGeneStateB(GeneBDE) = 2;
% Now create the expression data using functions
extGeneATotal = [extGeneA, extGeneC];

dataANonControl = CreateExpression(gGeneStateA, extGeneATotal, gRelatedGenes, ...
    dataAControl, dataANonControl, fixedSigma,  intGap, extGap);

extGeneBTotal = [extGeneB, extGeneC];
dataBNonControl = CreateExpression(gGeneStateB, extGeneBTotal, gRelatedGenes, ...
    dataBControl, dataBNonControl, fixedSigma,  intGap, extGap);


% Once the data is prepated you can save it.
save(outFile, 'dataAControl','dataANonControl', 'dataBControl', ...
    'dataBNonControl', 'intGeneA','extGeneA', ...
    'intGeneB', 'extGeneB', 'intGeneC', 'extGeneC', 'intGap', 'extGap');

end

function dataNonControl = CreateExpression(gGeneState, extGene, gRelatedGenes, ...
    dataControl, dataNonControl, fixedSigma, intGap, extGap)

% In the probArr the first position indicates extGeneA, extGeneB, extGeneC


% First normalize the non-control with respect to the control data points.
[lrow, lcol] = size(dataControl);

for j=1:lrow
	df = mean(dataNonControl(j,:)) - mean(dataControl(j,:));
	dataNonControl(j,:) = dataNonControl(j,:) - df;
end

genCount =0;

for i =1:size(gRelatedGenes)
    genCount = genCount +1;
    
    %display(['Data generated for geneCount: ', num2str(genCount)]);
    localGene = gRelatedGenes(i);
    
    localGeneState = gGeneState(localGene);
    
    if 2 == localGeneState
        
    	if 1 == sum(extGene == localGene)
        	toss = rand();
        	if toss < 0.5
            	 dataNonControl(localGene,:) = dataNonControl(localGene,:) - extGap * fixedSigma;
        	else
            	 dataNonControl(i,:) = dataNonControl(,:) + extGap * fixedSigma;
        	end
    	else
        	toss = rand();
            if toss < 0.5
               dataNonControl(i,:) = dataNonControl(i,:) - intGap * fixedSigma;
            else
               dataNonControl(i,:) = dataNonControl(i,:) + intGap * fixedSigma;
        	end
    	end
       	       
    end
        
end

end

function [intGeneA, extGeneA, intGeneB, extGeneB, intGeneC, extGeneC] = ...
    initializeGenes(relProg, maxNeigh, totalSeed, countArr, ...
    gRelatedGenes, gGeneState, inProb)
% start with a random node from the list

len = length(gRelatedGenes);
relatedSize = gRelatedGenes(len);

% Now do a simple breadth first search towards the righthand side as 
% the logic is that only the right hand side can be affected. 
% Try to find upto maxNeigh number of neighbors of startSeed that are
% not DE till now.
localSeed =0;
% extGene contains the external candidates.
intGeneA = [];
intGeneB = [];
intGeneC = [];

extGeneA = [];
extGeneB = [];
extGeneC = [];

ucount =0;
while localSeed < totalSeed
    % Start with a random one from the list of related genes
    startSeed = floor(1 + (relatedSize - 1)* rand());
	
    localRightC = relProg.getRightNeigh(startSeed);
    localCLen = length(localRightC);
    if 0 == localCLen 
        continue;
     end
        
    % Check if the seed is equally expressed.
    if 1 == gGeneState(startSeed)
        
        intGene = [];
        extGene = [];

        localSeed = localSeed + 1;
        % Make the gene differentially expressed.
        gGeneState(startSeed) = 2;
        % Also include the startSeed into the external structure.
        
        extGene = [extGene startSeed];
       
        thisCount =1;
        sibs =[];
        parent = startSeed;
        %display(['DE gene:', num2str(parent), ' Count:', num2str(localDE)]);
        loopCount =0;
        while thisCount < maxNeigh
            % Push the right hand children of parent into queue.
            rightChildren = relProg.getRightNeigh(parent);
            sibs = [sibs' rightChildren']';
            % Check if the queue is empty. If it is empty then break the
            % loop.
            if true == isempty(sibs)
                break;
            end
            
            % Get the next element from the queue. We are sure from the
            % last line that the queue is not empty.
            parent = sibs(1);
            sibs(1) = [];
            
            % If the parent is not DE then make it DE.
            if 1 == gGeneState(parent)
                % Do some tricks 
                % Check the number of its predecessors that are already DE.
                % Based on that make it DE with some probability.
                
                leftPred = relProg.getLeftNeigh(parent);
                % Calculate the number of DE predecessor
                leftDECount = sum(gGeneState(leftPred) == 2);
                
                localProb = 1 - (1 - inProb)^leftDECount;
                localRand = rand();
                
                if localRand > localProb
                    ucount = ucount + 1;
                    display(['Unsuccessful attempt:', num2str(parent),...
                        ' ucount: ', num2str(ucount)]);
                    continue;
                end
                
                intGene = [intGene parent];
                
                gGeneState(parent) = 2;
                %display(['DE gene:', num2str(parent), ' Count:', num2str(localDE)]);
                % increase the value of the iteration variables and check
                % if they have reached the limit.
                %localDE = localDE + 1;
                thisCount = thisCount + 1;
                if thisCount == maxNeigh
                    break;
                end
                
            else
                loopCount = loopCount +1;
                if loopCount == 10 
                    break
                end
            end
            
                  
        end
        
        % Get a type of dataset A, B or Common based on the probability
        % array. Do not make it random. It is already random. This ensures
        % that at least the promarily affected genes are not working that
        % way.
        if  localSeed <= countArr(1)
            % Assign it to intGeneComm and extGeneComm
            intGeneC = [intGeneC intGene];
            extGeneC = [extGeneC extGene];
            
        elseif localSeed <= countArr(2)
            % Assign it to intGeneA and extGeneA
            intGeneA = [intGeneA intGene];
            extGeneA = [extGeneA extGene];
            
        else
            % Assign it to intGeneB and exeGeneB
            intGeneB = [intGeneB intGene];
            extGeneB = [extGeneB extGene];
        end
                
    end

end
end

function rGenes = getRelatedGenes(localRelations)
% do a union and unique
rGenes = unique(union(localRelations(:,1), localRelations(:,2)));

end


function X = myrandn( patientNum )
X = randn(1, patientNum);

%Z = (X - mean(X))/std(X);
end

