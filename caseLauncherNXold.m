% Case launcher
clear; close all; 
close all


txt='mallas/comparacion.dat';qmax=1/25e-3^2;




msh=impq8(txt);
t=2e-3;

% qmax=100;

%% Preprocess

elementType='Q8';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
% loadCase='Uniform';         %'Uniform' 'Constant Bending' 'Variable Bending'

% Mesh generation
% [elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=[msh.elem.nod(:,1:4),msh.cord(1,[1 2]),1,1]%quadrilateralDomainMeshGenerator(elementType,'Straight',5,2,2,2,0,0);
elementNodesArray=msh.elem.nod;
nodesPositionArray=msh.cord(:,[1 2]);
% Deformed plot
meshPlot(elementNodesArray,nodesPositionArray,'b','No');

nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

Lwall=find(msh.cord(:,1)==min(msh.cord(:,1)));
Rwall=find(msh.cord(:,1)==max(msh.cord(:,1)));
Swall=find(msh.cord(:,2)==min(msh.cord(:,2)));

midline=find(msh.cord(:,1)==0);
% 
% Lwall=[84 211:222 447:448 1572:1580 1586:1588 1919:1920];
% Rwall=[687:695 701:703 816:818 1178:1187 1313 1364 1467:1468];
Lwall=sortrows([Lwall msh.cord(Lwall,:)],3);Lwall=Lwall(:,1);
Rwall=sortrows([Rwall msh.cord(Rwall,:)],3);Rwall=Rwall(:,1);
midline=sortrows([midline,msh.cord(midline,:)],3,'descend');midline=midline(:,1);


% SCLnod=midline(20:-1:1);
SCLnod=Rwall;


% M=@(y) -y*0.01;

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,200e9,0.3);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
% boundaryConditionsArray(vertexNodes(1),1:2) = true;
%boundaryConditionsArray(vertexNodes(3),1) = true;
% boundaryConditionsArray(sideNodes(end,find(sideNodes(end,:))),1) = true;

% Load definition
distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
pointLoadsArray(Lwall,1)=distributedLoad('uniforme',Lwall,-qmax,nodesPositionArray);
% pointLoadsArray(Rwall)=-M(msh.cord(Rwall,2));

%                 pointLoadsArray (sideNodes(2,1:size(find(sideNodes(2,:)),2)),1) = 1/3;
%                 pointLoadsArray (sideNodes(2,2:2:size(find(sideNodes(2,:)),2)),1) = 2/3;
%                 pointLoadsArray (vertexNodes(2),1) = 1/6;
%                 pointLoadsArray (vertexNodes(4),1) = 1/6;

nodiz=Lwall((msh.cord(Lwall,2)==0));
nodde=Rwall((msh.cord(Rwall,2)==0));
boundaryConditionsArray([nodde],[1 2])=true;
boundaryConditionsArray([Rwall],[1])=true;
boundaryConditionsArray([Swall],[1])=true;


%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,t);

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
[elementStressAtNodes]=stressRecext(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);

[r,c]=find(msh.elem.nod==SCLnod(end));
S=squeeze(elementStressAtNodes(r,c,:));
Svmmax=sqrt(S(1)^2+S(2)^2+3*S(3)^2);
Svm_b=zeros(length(SCLnod),1);
Svm_m=zeros(length(SCLnod),1);
Sfem=zeros(length(SCLnod),1);
% for i=1:length(SCLnod)
[Svm_m,Svm_b,Svm_F,Sfem,Sm,Sb,SF,Sfemx]=SCL(elementStressAtNodes,msh,SCLnod);
% end

%%
magnificationFactor=0.1;

% Deformed plot
% meshPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','Yes');

% Stresses plot
bandPlot(elementNodesArray,nodesPositionArray+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));