%% Case launcher
clear; close all; 
set(0,'DefaultFigureWindowStyle','docked')

E = 200e3; nu = 0.3;
pressureNormal = 100;

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension

[elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',1,1,10,10,0,0.05);
nodesPositionArray(:,1) = nodesPositionArray(:,1) + 1;

% Mesh plot
figure
meshPlot(elementNodesArray,nodesPositionArray,'b','Yes');
drawnow

nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
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

% Pressure
pointLoadsArray = distributedLoad_3(elementType,Lside,pointLoadsArray,nodesPositionArray,pressureNormal);
% pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodesPositionArray,-pressureNormal);

%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix);

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
[elementStressExtrapolated,elementStressAtGaussPoints]=stressRecext(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% [elementStressAtNodes]=stressRecovery_2(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=10;

% Deformed plot
figure
% meshPlot(elementNodesArray,nodesPositionArray,'b','No');
meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');

% Stresses plot
figure
meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,1)));
title('Tensiones en los nodos [Mpa]')

figure
meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressExtrapolated(:,:,1)));
title('Tensiones extrapoladas desde los ptos. de Gauss [Mpa]')

% Sxx = squeeze(elementStressAtNodes(:,:,1))'; Syy = squeeze(elementStressAtNodes(:,:,2))'; Szz = squeeze(elementStressAtNodes(:,:,3))'; Szx = squeeze(elementStressAtNodes(:,:,4))';
% vonMisesStress = 1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 );
% figure
% meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% bandPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',vonMisesStress);
% title('VM [Mpa]')

% Displacement Matrix
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm