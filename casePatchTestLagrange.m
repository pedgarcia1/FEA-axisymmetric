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

[elements1,nodes1,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,nElementsR,nElementsZ,0,distorsion);
vertexNodes2 = vertexNodes2 + size(elements1,1); 
sideNodes2 = sideNodes2 + size(nodes1,1);
nodes2(:,1) = nodes2(:,1) + b;
elements = [elements1;elements2+size(nodes1,1)];
nodes = [nodes1;nodes2];

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
% sideNodes1(2,:) - sideNodes2(4,:) = intereferencia
C = zeros(nConstraints,nTotalDof);
rightSideDOF = convertNode2Dof(sideNodes(2,:),nDimensions);
leftSideDOF = convertNode2Dof(sideNodes2(4,:),nDimensions);
n = 1;
for i = 1:2:length(rightSideDOF)
    C(n,rightSideDOF(i)) = 1*precond;
    C(n,leftSideDOF(i)) = -1*precond;
    n = n + 1;
end

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elements,nodes,constitutiveMatrix);

K = [stiffnessMatrix C'
    C zeros(nConstraints,nConstraints)];

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;
isFreeC = isFree; isFreeC((end+1):(end+nConstraints)) = true;


% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';
loadsVector((end+1):(end+nConstraints)) = 1*interferencia*precond; % interferencia????


% Equation solving
displacementsReducedVector = K(isFreeC,isFreeC)\loadsVector(isFreeC);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector(1:end-nConstraints);
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