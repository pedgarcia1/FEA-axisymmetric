%% Case launcher
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 103.1; c = b + 78.1; h = 600;
eEsf = 103.1;
interferencia = 0.31; precond = 1e7;
nElementsZ = 5; nElementsR = 5; distorsion = 0;   
planeStrainFlag = 0;

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

%% mallador
% cilindros zunchado
[elements1,nodes1,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,nElementsR,nElementsZ,0,distorsion);
nNodes1 = size(nodes1,1); nElements1 = size(elements1,1); nNodes2 = size(nodes2,1);
vertexNodes2 = vertexNodes2 + nNodes1; 
sideNodes2 = sideNodes2 + size(nodes1,1);
nodes2(:,1) = nodes2(:,1) + b - interferencia ;
msh.elements = elements1;
msh.nodes = nodes1;
msh.sideNodes = sideNodes;
msh.vertexNodes = vertexNodes1;

% esfera superior
meshAngle = 90;
meshInnerRadious = a;
meshLength = meshInnerRadious*meshAngle*pi/180;

[esf.elements,esf.nodes,esf.vertexNodes,esf.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,eEsf,5,nElementsR,meshAngle,0);

esf.nodes(:,2) = esf.nodes(:,2) + h;
meshPlot(esf.elements,esf.nodes,'r','No')

% merge mallas
tol = 0.001;
[msh] = mergeElAr(msh,esf,tol);

% Mesh plot
figure; meshPlot(msh.elements,msh.nodes,'b','Yes');


nodes = msh.nodes; elements = msh.elements;

nElements=size(elements,1);    %Number of elements
nNodes=size(nodes,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs
nConstraints = length(sideNodes(2,:));

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
if planeStrainFlag
    boundaryConditionsArray(sideNodes(1,:),2)=true;
    boundaryConditionsArray(sideNodes(3,:),2)=true;
    boundaryConditionsArray(sideNodes2(1,:),2)=true;
    boundaryConditionsArray(sideNodes2(3,:),2)=true;
else
%     boundaryConditionsArray(sideNodes(1,:),2)=true;
%     boundaryConditionsArray(sideNodes2(1,:),2)=true;
    nodosConZ0 = find(ismembertol(nodes(:,2),0));
    nodosConR0 = find(ismembertol(nodes(:,1),0));
    boundaryConditionsArray(nodosConR0,1) = false;
    boundaryConditionsArray(nodosConZ0,2) = false;
end

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);
pointLoadsArray = distributedLoad_3(elementType,Lside,pointLoadsArray,nodes,pressureNormal);
% pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodesPositionArray,-pressureNormal);

%% Solver

% Constraints
% - u1 + u2 = interferencia
% - sideNodes1(2,:) + sideNodes2(4,:) = intereferencia
C = zeros(nConstraints,nTotalDof);
rightSideDOF = convertNode2Dof(sideNodes(2,:),nDimensions);
leftSideDOF = convertNode2Dof(sideNodes2(4,:),nDimensions);
n = 1;
for i = 1:2:length(rightSideDOF)
    C(n,rightSideDOF(i)) = -1*precond; % u1
    C(n,leftSideDOF(i)) = 1*precond; % u2
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
loadsVector((end+1):(end+nConstraints)) = interferencia*precond; % interferencia

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
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm
