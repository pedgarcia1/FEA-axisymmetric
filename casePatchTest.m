%% Case launcher
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 1; b = 2;  h = 1;
nElementsZ = 2; nElementsR = 2;  distorsion = 0;   
numbering = 'Yes'; % Yes/No
tol = 0.1; 
caso = 2;

%% Preprocess

elementType = 'Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType = 'Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions = 2;              % Problem dimension

[elements,nodes,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes(:,1) = nodes(:,1) + a;
figure; hold on; meshPlot(elements,nodes,'b',numbering); xlim([0 b]);

nElements=size(elements,1);    %Number of elements
nNodes=size(nodes,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs
nConstraints = 0;

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
switch caso
    case 1
        boundaryConditionsArray(sideNodes(1,:),2) = true;
    case 2
        boundaryConditionsArray(sideNodes(1,:),2) = true;
        boundaryConditionsArray(sideNodes(3,:),2) = true;
    case 3
        boundaryConditionsArray(sideNodes(1,:),2) = true;
        boundaryConditionsArray(sideNodes(2,:),1) = true;        
        boundaryConditionsArray(sideNodes(4,:),1) = true; 
end

scatter(nodes(boundaryConditionsArray(:,1),1),nodes(boundaryConditionsArray(:,1),2),'filled','red');
scatter(nodes(boundaryConditionsArray(:,2),1),nodes(boundaryConditionsArray(:,2),2),'filled','green');

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);

switch caso
    case 1
        pointLoadsArray = distributedLoad_3(elementType,Lside,1,pointLoadsArray,nodes,pressureNormal);
        pointLoadsArray = distributedLoad_3(elementType,Rside,1,pointLoadsArray,nodes,-pressureNormal);
    case 2
        pointLoadsArray = distributedLoad_3(elementType,Lside,1,pointLoadsArray,nodes,pressureNormal);
        pointLoadsArray = distributedLoad_3(elementType,Rside,1,pointLoadsArray,nodes,-pressureNormal);
    case 3
        pointLoadsArray = distributedLoad_3(elementType,sideNodes(3,:),2,pointLoadsArray,nodes,-pressureNormal);
end
scatter(nodes(pointLoadsArray(:,1) ~= 0,1),nodes(pointLoadsArray(:,1) ~= 0,2),'filled','yellow')
% quiver(nodes(Lside,1),nodes(Lside,2),pointLoadsArray(Lside,1),pointLoadsArray(Lside,2))
% quiver(nodes(Rside,1),nodes(Rside,2),pointLoadsArray(Rside,1),pointLoadsArray(Rside,2))

%% Solver

% Constraints
% - u1 + u2 = interferencia
% - sideNodes1(2,:) + sideNodes2(4,:) = intereferencia
C = zeros(nConstraints,nTotalDof);
% rightSideDOF = convertNode2Dof(sideNodes(2,:),nDimensions);
% leftSideDOF = convertNode2Dof(sideNodes2(4,:),nDimensions);
% n = 1;
% for i = 1:2:length(rightSideDOF)
%     C(n,rightSideDOF(i)) = -1*precond; % u1
%     C(n,leftSideDOF(i)) = 1*precond; % u2
%     n = n + 1;
% end

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elements,nodes,constitutiveMatrix);
K = stiffnessMatrix;

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;
isFreeC = isFree; isFreeC((end+1):(end+nConstraints)) = true;

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = K(isFreeC,isFreeC)\loadsVector(isFreeC);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector(1:end-nConstraints);
% lagrangeMultipliers = displacementsReducedVector(end-nConstraints+1:end);

%% Postprocess
%Stress recovery
[elementStressExtrapolated,elementStressAtGaussPoints]=stressRecext(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
[elementStressAtNodes]=stressRecovery_2(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=30;
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm

%% Plots tension y malla deformada
% Deformed plot
figure; hold on; xlim([0 b]);
meshPlot(elements,nodes,'b','No');
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','r','No');

% Stresses plot
figure;
subplot(2,2,1); title('\sigma_{r}','Interpreter','tex')
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,1)),[-150 0]);
subplot(2,2,2); title('\sigma_{\theta}','Interpreter','tex')
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,2)),[-150 0]);
subplot(2,2,3); title('\sigma_{z}','Interpreter','tex')
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,3)),[-150 0]);
subplot(2,2,4); title('\tau_{zr}','Interpreter','tex')
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,4)),[-150 0]);



% % Von Mises Plot
% Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
% vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));
% figure; subplot(1,2,1); title('Int Svm [Mpa]')
% meshPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])','b','No');
% bandPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])',vonMisesStress(1:nNodes1,:));
%  subplot(1,2,2); title('Out Svm [Mpa]')
% meshPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])','b','No');
% bandPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])',vonMisesStress(nElements1+1:end,:));
% 
% % toda la malla VM Plot
% figure; title('Svm [MPa]')
% meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])','b','No');
% bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])',vonMisesStress);



