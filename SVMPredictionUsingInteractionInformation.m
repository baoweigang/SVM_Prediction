clear;
close all;
clc;
monitorPositions = get(0,'MonitorPositions');
set(0,'DefaultFigurePosition',monitorPositions(2,:));

saveOutputs             = false;
graphdpi                = 200;
maxStimsInSeq           = 7;
maxStimTime             = 1500; % 1.5 seconds per image, so 1500 ms
binTimeWidth            = 20;   % 20ms bins
binsInStimulus          = (maxStimTime/binTimeWidth);

colorLabels = cell(maxStimsInSeq,1);
colorLabels{1} = 'Varying LS';
colorLabels{2} = 'Varying NE';
colorLabels{3} = 'Varying TH';
colorLabels{4} = 'Object';
colorLabels{5} = 'Constant LS';
colorLabels{6} = 'Constant NE';
colorLabels{7} = 'Constant TH';

load('C:\Workspace_Rohit\Data Analysis\Outputs\SocialExperience\Neuron Coactivity Firing 20ms\Session01_Neuron_Coactivity_Patterns_20ms_6Cells.mat');
load('C:\Workspace_Rohit\Data Analysis\Outputs\SocialExperience\Neuron Coactivity Firing 20ms\stimulusBins.mat');
outputFolder ='C:\Workspace_Rohit\Data Analysis\Outputs\SocialExperience\StimulusStrings\';

nSessions = length(neuronCoactivity);

codewordArrays = cell(nSessions,1);
numberOfCells  = double(zeros(nSessions,1));
nUniqueWords   = double(zeros(nSessions,1));
for sessIdx = 1:nSessions
    codewordArrays{sessIdx} = neuronCoactivity(sessIdx).Codeword_Array;
    numberOfCells(sessIdx)  = log(size(neuronCoactivity(sessIdx).Image_Coactivity_Array,2))/log(2);    
    nUniqueWords(sessIdx)   = length(unique(codewordArrays{sessIdx}));
end

for sessIdx = 1:1%nSessions
    codewordArray = codewordArrays{sessIdx};
    nCells        = 6;%numberOfCells(sessIdx);       
    maxStates     = 2^nCells;
    maxBins       = length(codewordArray);
    
    tmpCounter = 0;    
    wordCounter = double(zeros(1000,maxStates,3));
    binStim1stSeq = double(zeros(maxBins,2));
    binStim2ndSeq = double(zeros(maxBins,2));
    binStim3rdSeq = double(zeros(maxBins,2));
    for binIdx = 2:maxBins
        if( (stimulusBins(1).Image_Color_Code(binIdx) == stimulusBins(1).Image_Color_Code(binIdx-1)) && (stimulusBins(1).Image_Color_Code(binIdx) > 0) )
            for wordIdx = 1:maxStates
                if(codewordArray(binIdx) == (wordIdx-1))
                    wordCounter(tmpCounter,wordIdx,stimulusBins(1).Image_Sequence_Code(binIdx)) = wordCounter(tmpCounter,wordIdx,stimulusBins(1).Image_Sequence_Code(binIdx)) + 1;
                    if(stimulusBins(1).Image_Sequence_Code(binIdx) == 1)
                        binStim1stSeq(binIdx,1) = stimulusBins(1).Image_Color_Code(binIdx);
                        binStim1stSeq(binIdx,2) = tmpCounter;
                    end
                    if(stimulusBins(1).Image_Sequence_Code(binIdx) == 2)
                        binStim2ndSeq(binIdx,1) = stimulusBins(1).Image_Color_Code(binIdx);
                        binStim2ndSeq(binIdx,2) = tmpCounter;
                    end
                    if(stimulusBins(1).Image_Sequence_Code(binIdx) == 3)
                        binStim3rdSeq(binIdx,1) = stimulusBins(1).Image_Color_Code(binIdx);
                        binStim3rdSeq(binIdx,2) = tmpCounter;
                    end
                end
            end
        elseif(stimulusBins(1).Image_Color_Code(binIdx) ~= stimulusBins(1).Image_Color_Code(binIdx-1))
            if(tmpCounter == 0)
                tmpCounter = tmpCounter + 1;
            else
                tmpCounter = tmpCounter + 0.5; % Yes, this seems weird, but there are two changes to the next stimulus (the falling edge of the current stimulus, and the rising edge of the next stimulus). So together this will increment the tmpCounter by 1, but the time the next stimulus rolls around.
            end
        end            
    end
    
    words1stSeq(1:1000,1:maxStates) = wordCounter(1:1000,1:maxStates,1);
    words2ndSeq(1:1000,1:maxStates) = wordCounter(1:1000,1:maxStates,2);
    words3rdSeq(1:1000,1:maxStates) = wordCounter(1:1000,1:maxStates,3);
    
    nStims1stSeq = nnz(sum(words1stSeq,2));
    nStims2ndSeq = nnz(sum(words2ndSeq,2));
    nStims3rdSeq = nnz(sum(words3rdSeq,2));
    
    words1stSeq = words1stSeq(1:nStims1stSeq,:);
    words2ndSeq = words2ndSeq(nStims1stSeq+1:nStims1stSeq+nStims2ndSeq,:);
    words3rdSeq = words3rdSeq(nStims1stSeq+nStims2ndSeq+1:nStims1stSeq+nStims2ndSeq+nStims3rdSeq,:);
    
    binStim1stSeqList = cell(nStims1stSeq,maxStimsInSeq);
    binStim2ndSeqList = cell(nStims2ndSeq,maxStimsInSeq);
    binStim3rdSeqList = cell(nStims3rdSeq,maxStimsInSeq);
    for binIdx = 1:maxBins
        if(binStim1stSeq(binIdx,1) > 0)            
            binStim1stSeqList{binStim1stSeq(binIdx,2),binStim1stSeq(binIdx,1)} = [binStim1stSeqList{binStim1stSeq(binIdx,2),binStim1stSeq(binIdx,1)},binIdx];
        end
        if(binStim2ndSeq(binIdx,1) > 0)            
            binStim2ndSeqList{(binStim2ndSeq(binIdx,2)-(nStims1stSeq)),binStim2ndSeq(binIdx,1)} = [binStim2ndSeqList{(binStim2ndSeq(binIdx,2)-(nStims1stSeq)),binStim2ndSeq(binIdx,1)},binIdx];
        end
        if(binStim3rdSeq(binIdx,1) > 0)            
            binStim3rdSeqList{(binStim3rdSeq(binIdx,2)-(nStims1stSeq+nStims2ndSeq)),binStim3rdSeq(binIdx,1)} = [binStim3rdSeqList{(binStim3rdSeq(binIdx,2)-(nStims1stSeq+nStims2ndSeq)),binStim3rdSeq(binIdx,1)},binIdx];
        end
    end
end

stimDist1st = double(zeros(maxStimsInSeq,1));
stimDist2nd = double(zeros(maxStimsInSeq,1));
stimDist3rd = double(zeros(maxStimsInSeq,1));
stimList1st = cell(maxStimsInSeq,1);
stimList2nd = cell(maxStimsInSeq,1);
stimList3rd = cell(maxStimsInSeq,1);
for stimIdx=1:maxStimsInSeq
    for stimNoIdx=1:nStims1stSeq    
        if(~isempty(binStim1stSeqList{stimNoIdx,stimIdx}))
            stimDist1st(stimIdx) = stimDist1st(stimIdx) + 1;
            stimList1st{stimIdx} = [stimList1st{stimIdx},stimNoIdx];
        end
    end
    for stimNoIdx=1:nStims2ndSeq
        if(~isempty(binStim2ndSeqList{stimNoIdx,stimIdx}))
            stimDist2nd(stimIdx) = stimDist2nd(stimIdx) + 1;            
            stimList2nd{stimIdx} = [stimList2nd{stimIdx},stimNoIdx];
        end
    end
    for stimNoIdx=1:nStims3rdSeq    
        if(~isempty(binStim3rdSeqList{stimNoIdx,stimIdx}))
            stimDist3rd(stimIdx) = stimDist3rd(stimIdx) + 1;
            stimList3rd{stimIdx} = [stimList3rd{stimIdx},stimNoIdx];
        end
    end
end

stimTest1st = cell(maxStimsInSeq,1);
stimTest2nd = cell(maxStimsInSeq,1);
stimTest3rd = cell(maxStimsInSeq,1);
for stimIdx = 1:maxStimsInSeq    
    % 1st Image Sequence
    clear tmpArray;
    tmpArray = stimList1st{stimIdx};
    nStimuli = length(tmpArray);
    if(nStimuli > 20)        
        nTestStims = 9; % Object so there's 3 times the stimulus
    else
        nTestStims = 3; 
    end
    
    clear tmpFlags;
    tmpFlags = false(nStimuli,1);
    tmpCounter = 0;
    while(tmpCounter < nTestStims)
        tmpNo = round((nStimuli-1)*rand(1) + 1);
        if(~tmpFlags(tmpNo))
            tmpFlags(tmpNo) = true;
            tmpCounter = tmpCounter + 1;
        end
    end
    stimTest1st{stimIdx} = tmpFlags;
    
    % 2nd Image Sequence
    clear tmpArray;
    tmpArray = stimList2nd{stimIdx};
    nStimuli = length(tmpArray);
    if(nStimuli > 20)        
        nTestStims = 9; % Object so there's 3 times the stimulus
    else
        nTestStims = 3; 
    end
    
    clear tmpFlags;
    tmpFlags = false(nStimuli,1);
    tmpCounter = 0;
    while(tmpCounter < nTestStims)
        tmpNo = round((nStimuli-1)*rand(1) + 1);
        if(~tmpFlags(tmpNo))
            tmpFlags(tmpNo) = true;
            tmpCounter = tmpCounter + 1;
        end
    end
    stimTest2nd{stimIdx} = tmpFlags;
    
    % 3rd Image Sequence
    clear tmpArray;
    tmpArray = stimList3rd{stimIdx};
    nStimuli = length(tmpArray);
    if(nStimuli > 20)        
        nTestStims = 9; % Object so there's 3 times the stimulus
    else
        nTestStims = 3; 
    end
    
    clear tmpFlags;
    tmpFlags = false(nStimuli,1);
    tmpCounter = 0;
    while(tmpCounter < nTestStims)
        tmpNo = round((nStimuli-1)*rand(1) + 1);
        if(~tmpFlags(tmpNo))
            tmpFlags(tmpNo) = true;
            tmpCounter = tmpCounter + 1;
        end
    end
    stimTest3rd{stimIdx} = tmpFlags;        
end

maxFeatures = 0; 
for cellIdx = 1:nCells
    maxFeatures = maxFeatures + nchoosek(nCells,cellIdx);
end    
maxFeatures = 2*maxFeatures + 2; % Entropy (1) and Information (2) = nchoosek; + Trace_Entropy (1) & Total_Correlation (2).

% Obtaining Mutual & Interaction Information For the 1st Image Sequence
featureVector1stSeq = double(NaN*ones(nStims1stSeq,maxFeatures));
for stimNoIdx = 1:nStims1stSeq
    clear tmpWordDistn;
    clear stimEntropy;
    tmpWordDistn(1:maxStates) = words1stSeq(stimNoIdx,1:maxStates);
    stimEntropy = interactionInformation(tmpWordDistn);%entropyPermutations_V2(tmpWordDistn);
    
    featCounter = 1;
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.All_Entropies{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector1stSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.Interaction_Information{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector1stSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    featureVector1stSeq(stimNoIdx,featCounter) = stimEntropy.Trace_Entropy;
    featCounter = featCounter + 1;
    featureVector1stSeq(stimNoIdx,featCounter) = stimEntropy.Total_Correlation;
    if(featCounter ~= maxFeatures)
        fprintf('Error: Something wrong, double check!\n');
    end
  
    fprintf('Computed Entropy & Information for stimNo %d of %d in the 1st Image Sequence.\n',stimNoIdx,nStims1stSeq);
end
classLabels1stSeq = double(zeros(nStims1stSeq,1));
for stimIdx = 1:maxStimsInSeq
    classLabels1stSeq(stimList1st{stimIdx}) = stimIdx;
end

testStimFlags = false(nStims1stSeq,1);
for stimIdx = 1:maxStimsInSeq
    clear tmpListLabels;
    clear tmpTestLabels;
    tmpListLabels = stimList1st{stimIdx};
    tmpTestLabels = stimTest1st{stimIdx};    
    testStimFlags(tmpListLabels(tmpTestLabels)) = true;    
end

trainFeatureVector1stSeq(:,1:maxFeatures) = featureVector1stSeq(testStimFlags==0,1:maxFeatures);
testtFeatureVector1stSeq(:,1:maxFeatures) = featureVector1stSeq(testStimFlags==1,1:maxFeatures);
trainClassLabels1stSeq(:,1)      = classLabels1stSeq(testStimFlags==0,1);
testtClassLabels1stSeq(:,1)      = classLabels1stSeq(testStimFlags==1,1);

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel1stSeq = fitcecoc(trainFeatureVector1stSeq,trainClassLabels1stSeq,'Learners',learnTemplate);
prediction1stSeq = predict(SVMModel1stSeq,testtFeatureVector1stSeq);

opPredict1stSeq = cell(maxStimsInSeq,1);
for tmpIdx=1:length(testtClassLabels1stSeq)
    fprintf('1st Image Sequence - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels1stSeq(tmpIdx),prediction1stSeq(tmpIdx));
    opPredict1stSeq{testtClassLabels1stSeq(tmpIdx)} = [opPredict1stSeq{testtClassLabels1stSeq(tmpIdx)},prediction1stSeq(tmpIdx)];
end

% Obtaining Mutual & Interaction Information For the 2nd Image Sequence
featureVector2ndSeq = double(NaN*ones(nStims2ndSeq,maxFeatures));
for stimNoIdx = 1:nStims2ndSeq
    clear tmpWordDistn;
    clear stimEntropy;
    tmpWordDistn(1:maxStates) = words2ndSeq(stimNoIdx,1:maxStates);
    stimEntropy = interactionInformation(tmpWordDistn);%entropyPermutations_V2(tmpWordDistn);
    
    featCounter = 1;
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.All_Entropies{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector2ndSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.Interaction_Information{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector2ndSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    featureVector2ndSeq(stimNoIdx,featCounter) = stimEntropy.Trace_Entropy;
    featCounter = featCounter + 1;
    featureVector2ndSeq(stimNoIdx,featCounter) = stimEntropy.Total_Correlation;
    if(featCounter ~= maxFeatures)
        fprintf('Error: Something wrong, double check!\n');
    end
 
    fprintf('Computed Entropy & Information for stimNo %d of %d in the 2nd Image Sequence.\n',stimNoIdx,nStims2ndSeq);
end
classLabels2ndSeq = double(zeros(nStims2ndSeq,1));
for stimIdx = 1:maxStimsInSeq
    classLabels2ndSeq(stimList2nd{stimIdx}) = stimIdx;
end

testStimFlags = false(nStims2ndSeq,1);
for stimIdx = 1:maxStimsInSeq
    clear tmpListLabels;
    clear tmpTestLabels;
    tmpListLabels = stimList2nd{stimIdx};
    tmpTestLabels = stimTest2nd{stimIdx};    
    testStimFlags(tmpListLabels(tmpTestLabels)) = true;    
end

trainFeatureVector2ndSeq(:,1:maxFeatures) = featureVector2ndSeq(testStimFlags==0,1:maxFeatures);
testtFeatureVector2ndSeq(:,1:maxFeatures) = featureVector2ndSeq(testStimFlags==1,1:maxFeatures);
trainClassLabels2ndSeq(:,1)      = classLabels2ndSeq(testStimFlags==0,1);
testtClassLabels2ndSeq(:,1)      = classLabels2ndSeq(testStimFlags==1,1);

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel2ndSeq = fitcecoc(trainFeatureVector2ndSeq,trainClassLabels2ndSeq,'Learners',learnTemplate);
prediction2ndSeq = predict(SVMModel2ndSeq,testtFeatureVector2ndSeq);

opPredict2ndSeq = cell(maxStimsInSeq,1);
for tmpIdx=1:length(testtClassLabels2ndSeq)
    fprintf('2nd Image Sequence - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels2ndSeq(tmpIdx),prediction2ndSeq(tmpIdx));
    opPredict2ndSeq{testtClassLabels2ndSeq(tmpIdx)} = [opPredict2ndSeq{testtClassLabels2ndSeq(tmpIdx)},prediction2ndSeq(tmpIdx)];
end

% Obtaining Mutual & Interaction Information For the 3rd Image Sequence
featureVector3rdSeq = double(NaN*ones(nStims3rdSeq,maxFeatures));
for stimNoIdx = 1:nStims3rdSeq
    clear tmpWordDistn;
    clear stimEntropy;
    tmpWordDistn(1:maxStates) = words3rdSeq(stimNoIdx,1:maxStates);
    
    stimEntropy = interactionInformation(tmpWordDistn);%entropyPermutations_V2(tmpWordDistn);
    
    featCounter = 1;
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.All_Entropies{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector3rdSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    for cellIdx = 1:nCells
        clear tmpArray;        
        tmpArray = stimEntropy.Interaction_Information{cellIdx};        
        for tmpIdx = 1:length(tmpArray)
            featureVector3rdSeq(stimNoIdx,featCounter) = tmpArray(tmpIdx);
            featCounter = featCounter + 1;
        end
    end
    featureVector3rdSeq(stimNoIdx,featCounter) = stimEntropy.Trace_Entropy;
    featCounter = featCounter + 1;
    featureVector3rdSeq(stimNoIdx,featCounter) = stimEntropy.Total_Correlation;
    if(featCounter ~= maxFeatures)
        fprintf('Error: Something wrong, double check!\n');
    end

    fprintf('Computed Entropy & Information for stimNo %d of %d in the 3rd Image Sequence.\n',stimNoIdx,nStims3rdSeq);
end
classLabels3rdSeq = double(zeros(nStims3rdSeq,1));
for stimIdx = 1:maxStimsInSeq
    classLabels3rdSeq(stimList3rd{stimIdx}) = stimIdx;
end

testStimFlags = false(nStims3rdSeq,1);
for stimIdx = 1:maxStimsInSeq
    clear tmpListLabels;
    clear tmpTestLabels;
    tmpListLabels = stimList3rd{stimIdx};
    tmpTestLabels = stimTest3rd{stimIdx};    
    testStimFlags(tmpListLabels(tmpTestLabels)) = true;    
end

trainFeatureVector3rdSeq(:,1:maxFeatures) = featureVector3rdSeq(testStimFlags==0,1:maxFeatures);
testtFeatureVector3rdSeq(:,1:maxFeatures) = featureVector3rdSeq(testStimFlags==1,1:maxFeatures);
trainClassLabels3rdSeq(:,1)      = classLabels3rdSeq(testStimFlags==0,1);
testtClassLabels3rdSeq(:,1)      = classLabels3rdSeq(testStimFlags==1,1);

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel3rdSeq = fitcecoc(trainFeatureVector3rdSeq,trainClassLabels3rdSeq,'Learners',learnTemplate);
prediction3rdSeq = predict(SVMModel3rdSeq,testtFeatureVector3rdSeq);

opPredict3rdSeq = cell(maxStimsInSeq,1);
for tmpIdx=1:length(testtClassLabels3rdSeq)
    fprintf('3rd Image Sequence - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels3rdSeq(tmpIdx),prediction3rdSeq(tmpIdx));
    opPredict3rdSeq{testtClassLabels3rdSeq(tmpIdx)} = [opPredict3rdSeq{testtClassLabels3rdSeq(tmpIdx)},prediction3rdSeq(tmpIdx)];
end

% By Identity
% 1st Image Sequence
trainClassLabels1stSeqByIdentity = double(zeros(length(trainClassLabels1stSeq),1));
trainClassLabels1stSeqByIdentity(trainClassLabels1stSeq < 4) = 1;
trainClassLabels1stSeqByIdentity(trainClassLabels1stSeq == 4) = 2;
trainClassLabels1stSeqByIdentity(trainClassLabels1stSeq > 4) = 3;

testtClassLabels1stSeqByIdentity = double(zeros(length(testtClassLabels1stSeq),1));
testtClassLabels1stSeqByIdentity(testtClassLabels1stSeq < 4) = 1;
testtClassLabels1stSeqByIdentity(testtClassLabels1stSeq == 4) = 2;
testtClassLabels1stSeqByIdentity(testtClassLabels1stSeq > 4) = 3;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel1stSeqByIdentity = fitcecoc(trainFeatureVector1stSeq,trainClassLabels1stSeqByIdentity,'Learners',learnTemplate);
prediction1stSeqByIdentity = predict(SVMModel1stSeqByIdentity,testtFeatureVector1stSeq);

opPredict1stSeqByIdentity = cell(3,1);
for tmpIdx=1:length(testtClassLabels1stSeqByIdentity)
    fprintf('1st Image Sequence By Identity - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels1stSeqByIdentity(tmpIdx),prediction1stSeqByIdentity(tmpIdx));
    opPredict1stSeqByIdentity{testtClassLabels1stSeqByIdentity(tmpIdx)} = [opPredict1stSeqByIdentity{testtClassLabels1stSeqByIdentity(tmpIdx)},prediction1stSeqByIdentity(tmpIdx)];
end

% 2nd Image Sequence
trainClassLabels2ndSeqByIdentity = double(zeros(length(trainClassLabels2ndSeq),1));
trainClassLabels2ndSeqByIdentity(trainClassLabels2ndSeq < 4) = 1;
trainClassLabels2ndSeqByIdentity(trainClassLabels2ndSeq == 4) = 2;
trainClassLabels2ndSeqByIdentity(trainClassLabels2ndSeq > 4) = 3;

testtClassLabels2ndSeqByIdentity = double(zeros(length(testtClassLabels2ndSeq),1));
testtClassLabels2ndSeqByIdentity(testtClassLabels2ndSeq < 4) = 1;
testtClassLabels2ndSeqByIdentity(testtClassLabels2ndSeq == 4) = 2;
testtClassLabels2ndSeqByIdentity(testtClassLabels2ndSeq > 4) = 3;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel2ndSeqByIdentity = fitcecoc(trainFeatureVector2ndSeq,trainClassLabels2ndSeqByIdentity,'Learners',learnTemplate);
prediction2ndSeqByIdentity = predict(SVMModel2ndSeqByIdentity,testtFeatureVector2ndSeq);

opPredict2ndSeqByIdentity = cell(3,1);
for tmpIdx=1:length(testtClassLabels2ndSeqByIdentity)
    fprintf('2nd Image Sequence By Identity - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels2ndSeqByIdentity(tmpIdx),prediction2ndSeqByIdentity(tmpIdx));
    opPredict2ndSeqByIdentity{testtClassLabels2ndSeqByIdentity(tmpIdx)} = [opPredict2ndSeqByIdentity{testtClassLabels2ndSeqByIdentity(tmpIdx)},prediction2ndSeqByIdentity(tmpIdx)];
end

% 3rd Image Sequence
trainClassLabels3rdSeqByIdentity = double(zeros(length(trainClassLabels3rdSeq),1));
trainClassLabels3rdSeqByIdentity(trainClassLabels3rdSeq < 4) = 1;
trainClassLabels3rdSeqByIdentity(trainClassLabels3rdSeq == 4) = 2;
trainClassLabels3rdSeqByIdentity(trainClassLabels3rdSeq > 4) = 3;

testtClassLabels3rdSeqByIdentity = double(zeros(length(testtClassLabels3rdSeq),1));
testtClassLabels3rdSeqByIdentity(testtClassLabels3rdSeq < 4) = 1;
testtClassLabels3rdSeqByIdentity(testtClassLabels3rdSeq == 4) = 2;
testtClassLabels3rdSeqByIdentity(testtClassLabels3rdSeq > 4) = 3;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel3rdSeqByIdentity = fitcecoc(trainFeatureVector3rdSeq,trainClassLabels3rdSeqByIdentity,'Learners',learnTemplate);
prediction3rdSeqByIdentity = predict(SVMModel3rdSeqByIdentity,testtFeatureVector3rdSeq);

opPredict3rdSeqByIdentity = cell(3,1);
for tmpIdx=1:length(testtClassLabels3rdSeqByIdentity)
    fprintf('3rd Image Sequence By Identity - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels3rdSeqByIdentity(tmpIdx),prediction3rdSeqByIdentity(tmpIdx));
    opPredict3rdSeqByIdentity{testtClassLabels3rdSeqByIdentity(tmpIdx)} = [opPredict3rdSeqByIdentity{testtClassLabels3rdSeqByIdentity(tmpIdx)},prediction3rdSeqByIdentity(tmpIdx)];
end

% By Stimulus
% 1st Image Sequence
trainClassLabels1stSeqByStimulus = double(zeros(length(trainClassLabels1stSeq),1));
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 1) = 1;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 5) = 1;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 2) = 2;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 6) = 2;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 3) = 3;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 7) = 3;
trainClassLabels1stSeqByStimulus(trainClassLabels1stSeq == 4) = 4;

testtClassLabels1stSeqByStimulus = double(zeros(length(testtClassLabels1stSeq),1));
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 1) = 1;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 5) = 1;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 2) = 2;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 6) = 2;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 3) = 3;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 7) = 3;
testtClassLabels1stSeqByStimulus(testtClassLabels1stSeq == 4) = 4;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel1stSeqByStimulus = fitcecoc(trainFeatureVector1stSeq,trainClassLabels1stSeqByStimulus,'Learners',learnTemplate);
prediction1stSeqByStimulus = predict(SVMModel1stSeqByStimulus,testtFeatureVector1stSeq);

opPredict1stSeqByStimulus = cell(4,1);
for tmpIdx=1:length(testtClassLabels1stSeqByStimulus)
    fprintf('1st Image Sequence By Stimulus - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels1stSeqByStimulus(tmpIdx),prediction1stSeqByStimulus(tmpIdx));
    opPredict1stSeqByStimulus{testtClassLabels1stSeqByStimulus(tmpIdx)} = [opPredict1stSeqByStimulus{testtClassLabels1stSeqByStimulus(tmpIdx)},prediction1stSeqByStimulus(tmpIdx)];
end

% 2nd Image Sequence
trainClassLabels2ndSeqByStimulus = double(zeros(length(trainClassLabels2ndSeq),1));
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 1) = 1;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 5) = 1;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 2) = 2;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 6) = 2;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 3) = 3;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 7) = 3;
trainClassLabels2ndSeqByStimulus(trainClassLabels2ndSeq == 4) = 4;

testtClassLabels2ndSeqByStimulus = double(zeros(length(testtClassLabels2ndSeq),1));
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 1) = 1;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 5) = 1;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 2) = 2;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 6) = 2;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 3) = 3;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 7) = 3;
testtClassLabels2ndSeqByStimulus(testtClassLabels2ndSeq == 4) = 4;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel2ndSeqByStimulus = fitcecoc(trainFeatureVector2ndSeq,trainClassLabels2ndSeqByStimulus,'Learners',learnTemplate);
prediction2ndSeqByStimulus = predict(SVMModel2ndSeqByStimulus,testtFeatureVector2ndSeq);

opPredict2ndSeqByStimulus = cell(4,1);
for tmpIdx=1:length(testtClassLabels2ndSeqByStimulus)
    fprintf('2nd Image Sequence By Stimulus - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels2ndSeqByStimulus(tmpIdx),prediction2ndSeqByStimulus(tmpIdx));
    opPredict2ndSeqByStimulus{testtClassLabels2ndSeqByStimulus(tmpIdx)} = [opPredict2ndSeqByStimulus{testtClassLabels2ndSeqByStimulus(tmpIdx)},prediction2ndSeqByStimulus(tmpIdx)];
end

% 3rd Image Sequence
trainClassLabels3rdSeqByStimulus = double(zeros(length(trainClassLabels3rdSeq),1));
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 1) = 1;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 5) = 1;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 2) = 2;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 6) = 2;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 3) = 3;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 7) = 3;
trainClassLabels3rdSeqByStimulus(trainClassLabels3rdSeq == 4) = 4;

testtClassLabels3rdSeqByStimulus = double(zeros(length(testtClassLabels3rdSeq),1));
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 1) = 1;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 5) = 1;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 2) = 2;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 6) = 2;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 3) = 3;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 7) = 3;
testtClassLabels3rdSeqByStimulus(testtClassLabels3rdSeq == 4) = 4;

% SVM Classification
learnTemplate = templateSVM('Standardize',1,'SaveSupportVectors',true);
SVMModel3rdSeqByStimulus = fitcecoc(trainFeatureVector3rdSeq,trainClassLabels3rdSeqByStimulus,'Learners',learnTemplate);
prediction3rdSeqByStimulus = predict(SVMModel3rdSeqByStimulus,testtFeatureVector3rdSeq);

opPredict3rdSeqByStimulus = cell(4,1);
for tmpIdx=1:length(testtClassLabels3rdSeqByStimulus)
    fprintf('3rd Image Sequence By Stimulus - Test Case %d: Actual = %d, Prediction = %d.\n',tmpIdx,testtClassLabels3rdSeqByStimulus(tmpIdx),prediction3rdSeqByStimulus(tmpIdx));
    opPredict3rdSeqByStimulus{testtClassLabels3rdSeqByStimulus(tmpIdx)} = [opPredict3rdSeqByStimulus{testtClassLabels3rdSeqByStimulus(tmpIdx)},prediction3rdSeqByStimulus(tmpIdx)];
end

% Displays
clc;
disp1stSeq = double(NaN*ones(maxStimsInSeq,9));
disp2ndSeq = double(NaN*ones(maxStimsInSeq,9));
disp3rdSeq = double(NaN*ones(maxStimsInSeq,9));
accuracy1stSeq = double(zeros(maxStimsInSeq,1));
accuracy2ndSeq = double(zeros(maxStimsInSeq,1));
accuracy3rdSeq = double(zeros(maxStimsInSeq,1));
for stimIdx = 1:maxStimsInSeq
    if(stimIdx ~= 4)
        disp1stSeq(stimIdx,1:3) = opPredict1stSeq{stimIdx};
        disp2ndSeq(stimIdx,1:3) = opPredict2ndSeq{stimIdx};
        disp3rdSeq(stimIdx,1:3) = opPredict3rdSeq{stimIdx};
    else
        disp1stSeq(stimIdx,1:9) = opPredict1stSeq{stimIdx};
        disp2ndSeq(stimIdx,1:9) = opPredict2ndSeq{stimIdx};
        disp3rdSeq(stimIdx,1:9) = opPredict3rdSeq{stimIdx};
    end 
    accuracy1stSeq(stimIdx) = length(find(disp1stSeq(stimIdx,:) == stimIdx))/sum(~isnan(disp1stSeq(stimIdx,:)));
    accuracy2ndSeq(stimIdx) = length(find(disp2ndSeq(stimIdx,:) == stimIdx))/sum(~isnan(disp2ndSeq(stimIdx,:)));
    accuracy3rdSeq(stimIdx) = length(find(disp3rdSeq(stimIdx,:) == stimIdx))/sum(~isnan(disp3rdSeq(stimIdx,:)));
end

disp1stSeqByIdentity = double(NaN*ones(3,9));
disp2ndSeqByIdentity = double(NaN*ones(3,9));
disp3rdSeqByIdentity = double(NaN*ones(3,9));
accuracy1stSeqByIdentity = double(zeros(3,1));
accuracy2ndSeqByIdentity = double(zeros(3,1));
accuracy3rdSeqByIdentity = double(zeros(3,1));
for stimIdx = 1:3    
    disp1stSeqByIdentity(stimIdx,1:9) = opPredict1stSeqByIdentity{stimIdx};
    disp2ndSeqByIdentity(stimIdx,1:9) = opPredict2ndSeqByIdentity{stimIdx};
    disp3rdSeqByIdentity(stimIdx,1:9) = opPredict3rdSeqByIdentity{stimIdx};
    
    accuracy1stSeqByIdentity(stimIdx) = length(find(disp1stSeqByIdentity(stimIdx,:) == stimIdx))/sum(~isnan(disp1stSeqByIdentity(stimIdx,:)));
    accuracy2ndSeqByIdentity(stimIdx) = length(find(disp2ndSeqByIdentity(stimIdx,:) == stimIdx))/sum(~isnan(disp2ndSeqByIdentity(stimIdx,:)));
    accuracy3rdSeqByIdentity(stimIdx) = length(find(disp3rdSeqByIdentity(stimIdx,:) == stimIdx))/sum(~isnan(disp3rdSeqByIdentity(stimIdx,:)));
end

disp1stSeqByStimulus = double(NaN*ones(4,9));
disp2ndSeqByStimulus = double(NaN*ones(4,9));
disp3rdSeqByStimulus = double(NaN*ones(4,9));
accuracy1stSeqByStimulus = double(zeros(4,1));
accuracy2ndSeqByStimulus = double(zeros(4,1));
accuracy3rdSeqByStimulus = double(zeros(4,1));
for stimIdx = 1:4
    if(stimIdx ~= 4)
        disp1stSeqByStimulus(stimIdx,1:6) = opPredict1stSeqByStimulus{stimIdx};
        disp2ndSeqByStimulus(stimIdx,1:6) = opPredict2ndSeqByStimulus{stimIdx};
        disp3rdSeqByStimulus(stimIdx,1:6) = opPredict3rdSeqByStimulus{stimIdx};
    else
        disp1stSeqByStimulus(stimIdx,1:9) = opPredict1stSeqByStimulus{stimIdx};
        disp2ndSeqByStimulus(stimIdx,1:9) = opPredict2ndSeqByStimulus{stimIdx};
        disp3rdSeqByStimulus(stimIdx,1:9) = opPredict3rdSeqByStimulus{stimIdx};
    end
    
    accuracy1stSeqByStimulus(stimIdx) = length(find(disp1stSeqByStimulus(stimIdx,:) == stimIdx))/sum(~isnan(disp1stSeqByStimulus(stimIdx,:)));
    accuracy2ndSeqByStimulus(stimIdx) = length(find(disp2ndSeqByStimulus(stimIdx,:) == stimIdx))/sum(~isnan(disp2ndSeqByStimulus(stimIdx,:)));
    accuracy3rdSeqByStimulus(stimIdx) = length(find(disp3rdSeqByStimulus(stimIdx,:) == stimIdx))/sum(~isnan(disp3rdSeqByStimulus(stimIdx,:)));
end

fprintf('1st Image Sequence\n');
disp1stSeq
accuracy1stSeq
fprintf('2nd Image Sequence\n');
disp2ndSeq
accuracy2ndSeq
fprintf('3rd Image Sequence\n');
disp3rdSeq
accuracy3rdSeq

fprintf('1st Image Sequence By Identity\n');
disp1stSeqByIdentity
accuracy1stSeqByIdentity
fprintf('2nd Image Sequence By Identity\n');
disp2ndSeqByIdentity
accuracy2ndSeqByIdentity
fprintf('3rd Image Sequence By Identity\n');
disp3rdSeqByIdentity
accuracy3rdSeqByIdentity

fprintf('1st Image Sequence By Stimulus\n');
disp1stSeqByStimulus
accuracy1stSeqByStimulus
fprintf('2nd Image Sequence By Stimulus\n');
disp2ndSeqByStimulus
accuracy2ndSeqByStimulus
fprintf('3rd Image Sequence By Stimulus\n');
disp3rdSeqByStimulus
accuracy3rdSeqByStimulus
