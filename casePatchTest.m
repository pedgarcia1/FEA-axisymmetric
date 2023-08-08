%% Case launcher
clear; close all; 
set(0,'DefaultFigureWindowStyle','docked')

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 1; b = 2; c = 3; h = 1;
interferencia = 0.4;
nElementsZ = 10; nElementsR = 10; distorsion = 0;

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension

[elements1,nodes1,vertexNodes1,sideNodes1]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,nElementsR,nElementsZ,0,distorsion);
vertexNodes2 = vertexNodes2 + size(elements1,1); 
sideNodes2 = sideNodes2 + size(elements1,1);
nodes2(:,1) = nodes2(:,1) + b - interferencia;
elements = [elements1;elements2];
nodes = [nodes1;nodes2];

% Mesh plot
figure; meshPlot(elements2,nodes2,'b','Yes'); drawnow;

nElements=size(elements,1);    %Number of elements
nNodes=size(nodes,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(sideNodes(1,:),2)=true;
boundaryConditionsArray(sideNodes(3,:),2)=true;

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);
pointLoadsArray = distributedLoad_3(elementType,Lside,pointLoadsArray,nodes,pressureNormal);
% pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodesPositionArray,-pressureNormal);

%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elements,nodes,constitutiveMatrix);

% Constraints
% sideNodes1(2,:) - sideNodes2(4,:) = intereferencia

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = stiffnessMatrix(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;

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