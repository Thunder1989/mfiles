clear
clc

% Make rain labels y, and binary month features X
% data = csvread('');
% y = data(:,1);
% data = data(:,2:end);

%fake data
data = randi(5,600,1000);
y = mod(1:size(data,1),10);
nInstances = 1;
nNodes = size(data,1);

%% Make edgeStruct
IDX = knnsearch(data,data,'k',4,'distance','euclidean');
nStates = length(unique(y));
adj = zeros(nNodes);
for i = 1:nNodes
    nn = IDX(i,:);
    for j = 1:length(nn)
        if nn(j) ~= i
            adj(i,nn(j)) = 1;
        end
    end    
end
% adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
nEdges = edgeStruct.nEdges;
maxState = max(nStates);
fprintf('edges constructed(paused)\n');
pause

%% Training (with node features, but no edge features)

% Make simple bias features
Xedge = ones(nInstances,1,nEdges);

% Make simple bias features
nFeatures = 12;
Xnode = zeros(nInstances,nFeatures,nNodes);
for m = 1:nFeatures
    Xnode(:,m,:) = 1;
end
Xnode = [ones(nInstances,1,nNodes) Xnode];
nNodeFeatures = size(Xnode,2);

% Make nodeMap
nodeMap = zeros(nNodes,maxState,nNodeFeatures,'int32');
for f = 1:nNodeFeatures
    nodeMap(:,1,f) = f;
end

% Make edgeMap
edgeMap = zeros(maxState,maxState,nEdges,'int32');
edgeMap(1,1,:) = nNodeFeatures+1;
edgeMap(2,1,:) = nNodeFeatures+2;
edgeMap(1,2,:) = nNodeFeatures+3;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Optimize
w = minFunc(@UGM_CRF_NLL,w,[],Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP)
fprintf('(paused)\n');
pause

%% Training (with edge features)

% Make edge features
sharedFeatures = 1:13;
Xedge = UGM_makeEdgeFeatures(Xnode,edgeStruct.edgeEnds,sharedFeattures);
nEdgeFeatures = size(Xedge,2);

% Make edgeMap
edgeMap = zeros(maxState,maxState,nEdges,nEdgeFeatures,'int32');
for edgeFeat = 1:nEdgeFeatures
    for s1 = 1:2
        for s2 = 1:2
            f = f+1;
            edgeMap(s1,s2,:,edgeFeat) = f;
        end
    end
end

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Optimize
UGM_CRF_NLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP);
w = minFunc(@UGM_CRF_NLL,w,[],Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP)
fprintf('(paused)\n');
pause

%% Do decoding/infence/sampling in learned model (given features)

% We will look at a case in December
i = 11;
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);

decode = UGM_Decode_Chain(nodePot,edgePot,edgeStruct)

[nodeBel,edgeBel,logZ] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure(1);
imagesc(samples')
title('Samples from CRF model (for December)');
fprintf('(paused)\n');
pause

%% Do conditional decoding/inference/sampling in learned model (given features)

clamped = zeros(nNodes,1);
clamped(1:2) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_Chain)
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_Chain)
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Chain);

figure(2);
imagesc(condSamples')
title('Conditional samples from CRF model (for December)');
fprintf('(paused)\n');
pause

%% Now see what samples in July look like

XtestNode = [1 0 0 0 0 0 0 1 0 0 0 0 0]; % Turn on bias and indicator variable for July
XtestNode = repmat(XtestNode,[1 1 nNodes]);
XtestEdge = UGM_makeEdgeFeatures(XtestNode,edgeStruct.edgeEnds,sharedFeatures);

[nodePot,edgePot] = UGM_CRF_makePotentials(w,XtestNode,XtestEdge,nodeMap,edgeMap,edgeStruct);

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure(3);
imagesc(samples')
title('Samples from CRF model (for July)');
fprintf('(paused)\n');
pause

%% Training with L2-regularization

% Set up regularization parameters
lambda = 10*ones(size(w));
lambda(1) = 0; % Don't penalize node bias variable
lambda(14:17) = 0; % Don't penalize edge bias variable
regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain);

% Optimize
w = zeros(nParams,1);
w = minFunc(regFunObj,w);
NLL = UGM_CRF_NLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)
