%% Case launcher
clear; close all; 
set(0,'DefaultFigureWindowStyle','docked')

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 69.1; c = b + 100.1; h = 500;
interferencia = 0.7427; precond = 1e7;
nElementsZ = 20; nElementsR = 20; distorsion = 0;

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

[elements,nodes,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes(:,1) = nodes(:,1) + a;

% Mesh plot
figure; meshPlot(elements,nodes,'b','No');

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


% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);
pointLoadsArray = distributedLoad_3(elementType,Lside,pointLoadsArray,nodes,pressureNormal);
pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodes,-51.8);

%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elements,nodes,constitutiveMatrix);

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
% lagrangeMultipliers = displacementsReducedVector(end-nConstraints+1:end);

%% Postprocess
%Stress recovery
[elementStressExtrapolated,elementStressAtGaussPoints]=stressRecext(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% [elementStressAtNodes]=stressRecovery_2(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=1;

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

Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));
figure
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',vonMisesStress);
title('VM [Mpa]')

% Displacement Matrix
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm