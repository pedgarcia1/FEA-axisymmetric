%% Case launcher
clear; close all; 
set(0,'DefaultFigureWindowStyle','docked')

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 101.1; c = b + 74.1; h = 250;
interferencia = 50;
nElementsZ = 20; nElementsR = 20; distorsion = 0;

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

[elements1,nodes1,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,nElementsR,nElementsZ,0,distorsion);
vertexNodes2 = vertexNodes2 + size(elements1,1); 
sideNodes2 = sideNodes2 + size(nodes1,1);
nodes2(:,1) = nodes2(:,1) + b - interferencia;
elements = [elements1;elements2+size(nodes1,1)];
nodes = [nodes1;nodes2];

% Mesh plot
figure; meshPlot(elements,nodes,'b','Yes');

nElements=size(elements,1);    %Number of elements
nNodes=size(nodes,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs
nConstraints = length(sideNodes(2,:));

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(sideNodes(1,:),2)=true;
boundaryConditionsArray(sideNodes(3,:),2)=true;
boundaryConditionsArray(sideNodes2(1,:),2)=true;
boundaryConditionsArray(sideNodes2(3,:),2)=true;

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);
pointLoadsArray = distributedLoad_3(elementType,Lside,pointLoadsArray,nodes,pressureNormal);
% pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodesPositionArray,-pressureNormal);

%% Solver

% Constraints
rightSideDOF = convertNode2Dof(sideNodes(2,:),nDimensions);
leftSideDOF = convertNode2Dof(sideNodes2(4,:),nDimensions);
leftSideDOF = leftSideDOF(1:2:end);
T = eye(nTotalDof-nConstraints,nTotalDof-nConstraints);
% n = 1;
% for i = rightSideDOF(1:2:end)-1
%     aux = T(1:i,:);
%     row = zeros(1,nTotalDof-nConstraints);
%     row(1,leftSideDOF(n)) = 1 + interferencia;
%     n = n + 1;
%     last = T(i+1:end,:);
%     T = [aux;row;last];
%     iPrev = i;
% end
n = 1;
leftSideDOF = [19 rightSideDOF(2:2:end)];
for i = rightSideDOF(1:2:end)-1
    aux = T(1:i,:);
    row = zeros(1,nTotalDof-nConstraints);
    row(1,leftSideDOF(n)) = 1 ;
    n = n + 1;
    last = T(i+1:end,:);
    T = [aux;row;last];
    iPrev = i;
end
isCondensed = ismember(1:nTotalDof,rightSideDOF(1:2:end));
isReleased = ~isCondensed;

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elements,nodes,constitutiveMatrix);

Kr = T'*stiffnessMatrix*T;

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';
Rr = T'*loadsVector;

% Equation solving
displacementsReducedVector = Kr(isFree(isReleased),isFree(isReleased))\Rr(isFree(isReleased));

% Reconstruction
Dr = zeros(nTotalDof-nConstraints,1);
Dr(isFree(isReleased)) = Dr(isFree(isReleased)) + displacementsReducedVector;
% lagrangeMultipliers = displacementsReducedVector(end-nConstraints+1:end);
displacementsVector = T*Dr;

%% Postprocess
%Stress recovery
[elementStressExtrapolated,elementStressAtGaussPoints]=stressRecext(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% [elementStressAtNodes]=stressRecovery_2(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=10;

% Deformed plot
figure
% meshPlot(elementNodesArray,nodesPositionArray,'b','No');
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');

% Stresses plot
figure
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,1)));
title('Tensiones en los nodos [Mpa]')

figure
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressExtrapolated(:,:,1)));
title('Tensiones extrapoladas desde los ptos. de Gauss [Mpa]')

% Sxx = squeeze(elementStressAtNodes(:,:,1))'; Syy = squeeze(elementStressAtNodes(:,:,2))'; Szz = squeeze(elementStressAtNodes(:,:,3))'; Szx = squeeze(elementStressAtNodes(:,:,4))';
% vonMisesStress = 1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 );
% figure
% meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% bandPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',vonMisesStress);
% title('VM [Mpa]')

% Displacement Matrix
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm