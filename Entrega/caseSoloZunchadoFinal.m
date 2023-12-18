%% Case launcher
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 69.1; c = b + 100.1; h = 500;
% eEsf = b-a; nElementsInRadius = 50; centroEsfera = [0 h];
interferencia = 0.1546; precond = 1e7;
nElementsZ = 59; nElementsR = 12;  distorsion = 0;   
planeStrainFlag = 0; numbering = 'No'; % Yes/No
tol = 0.1; 

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

%% mallador
[elements1,nodes1,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,8,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,12,nElementsZ,0,distorsion);
nNodes1 = size(nodes1,1); nElements1 = size(elements1,1); nNodes2 = size(nodes2,1);
vertexNodes2 = vertexNodes2 + nNodes1; 
sideNodes2 = sideNodes2 + size(nodes1,1);
nodes2(:,1) = nodes2(:,1) + b - interferencia ;
msh.elements = [elements1;elements2+size(nodes1,1)];
elements = msh.elements;
msh.nodes = [nodes1;nodes2];
nodes = msh.nodes;
msh.sideNodes = sideNodes;
msh.vertexNodes = vertexNodes1;
%% 

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
    boundaryConditionsArray(nodosConR0,1) = true;
    boundaryConditionsArray(nodosConZ0,2) = true;
end

figure; hold on; meshPlot(elements,nodes,'b',numbering);
scatter(nodes(boundaryConditionsArray(:,1),1),nodes(boundaryConditionsArray(:,1),2),'filled','red');
scatter(nodes(boundaryConditionsArray(:,2),1),nodes(boundaryConditionsArray(:,2),2),'filled','green');

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
Rside = sideNodes(2,:);
Lside = sideNodes(4,:);
pointLoadsArray = distributedLoad_3(elementType,Lside,1,pointLoadsArray,nodes,pressureNormal);
% pointLoadsArray = distributedLoad_3(elementType,Rside,pointLoadsArray,nodesPositionArray,-pressureNormal);

scatter(nodes(pointLoadsArray(:,1) ~= 0 | pointLoadsArray(:,2) ~= 0,1),nodes(pointLoadsArray(:,1) ~= 0 | pointLoadsArray(:,2) ~= 0,2),'filled','yellow')

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
[elementStressAtNodes]=stressRecovery_2(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=1;
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm
%% Solucion teorica zunchado
if planeStrainFlag
    E = E/(1-nu^2);
    nu = nu/(1-nu);
    fprintf("Caso plane strain \n")
else
    fprintf("Caso plane stress \n")
end

stressNX = [];
stressNX = importResultadosNX('soloZunchadoStress.csv');

% pInterferencia = abs(mean(mean(elementStressAtGaussPoints([90 301],:,1))))
outCilElements = find(any(ismember(elements,sideNodes2(4,:))'));
intCilElements = find(any(ismember(elements,sideNodes(2,:))'));
pInterferenciaInt = abs(mean(mean(elementStressAtGaussPoints(intCilElements,:,1))))
pInterferenciaOut = abs(mean(mean(elementStressAtGaussPoints(outCilElements,:,1))))
pIntAvg = mean([pInterferenciaOut pInterferenciaInt])
pInterferencia = 51.8;

C1= @(a,b,pInt,pOut) (a^2*pInt-b^2*pOut)/(b^2-a^2);
C2= @(a,b,pInt,pOut) (pInt-pOut)*a^2*b^2/(b^2-a^2);

sigmaR1=@(r)    C1(a,b,pressureNormal,pInterferencia)-C2(a,b,pressureNormal,pInterferencia)./r.^2;
sigmaTita1=@(r) C1(a,b,pressureNormal,pInterferencia)+C2(a,b,pressureNormal,pInterferencia)./r.^2;

sigmaR2=@(r)    C1(b,c,pInterferencia,0)-C2(b,c,pInterferencia,0)./r.^2;
sigmaTita2=@(r) C1(b,c,pInterferencia,0)+C2(b,c,pInterferencia,0)./r.^2;

figure; subplot(1,2,1); hold on; title('\sigma_{r}','Interpreter','tex','FontSize',16); grid
zPlot = 500;
log = find(abs(nodes(:,2) - zPlot) < 10);
for i = 1:size(log,1)
    eleLog = ismember(elements,log(i));
    iEles = find(any(eleLog')');
    for n = 1:size(iEles)
        auxStress(n) = elementStressExtrapolated(iEles(n),eleLog(iEles(n),:)',1);
    end
    h(1) = scatter(nodes(log(i),1),mean(auxStress));

    auxStress = [];
end
    resNXRR = importdata('resultadosZunchadoStressRR.txt');
    h(2) = scatter(resNXRR(:,2)+a,resNXRR(:,3),'x');
aux = plot(a:0.1:b,sigmaR1(a:0.1:b),'r',b:0.1:c,sigmaR2(b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','best','FontSize',16)
xlabel('r [mm]','FontSize',16)
ylabel('MPa','FontSize',16)

subplot(1,2,2); hold on; title('\sigma_{\theta}','Interpreter','tex','FontSize',16); grid
log = find(abs(nodes(:,2) - zPlot) < 10);
for i = 1:size(log,1)
    eleLog = ismember(elements,log(i));
    iEles = find(any(eleLog')');
    for n = 1:size(iEles)
        auxStress(n) = elementStressExtrapolated(iEles(n),eleLog(iEles(n),:)',2);
    end
    h(1) = scatter(nodes(log(i),1),mean(auxStress));
    
    auxStress = [];
end

resNXTT = importdata('resultadosZunchadoStressTT.txt');
    h(2) = scatter(resNXTT(:,2)+a,resNXTT(:,3),'x');


aux = plot(a:0.1:b,sigmaTita1(a:0.1:b),'r',b:0.1:c,sigmaTita2(b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','best','FontSize',16)
xlabel('r [mm]','FontSize',16)
ylabel('MPa','FontSize',16)
%%
uTeorico= @(a,b,pInt,pOut,r) ((1-nu)/E).* C1(a,b,pInt,pOut) .* r + ((1+nu)/E).*  C2(a,b,pInt,pOut) ./r;
uNX = importResultadosNXZunchado('soloZunchadoU.csv');

figure; hold on; title('u','FontSize',16,'Interpreter','latex'); grid
log = find(abs(nodes(:,2) - zPlot) < 10);
aux = scatter(nodes(log,1),displacementsMatrix(log,1));
h(1) = aux(1);
aux = scatter(uNX.RCoord,uNX.R,'x');
h(2) = aux(1);
aux = plot(a:0.1:b,uTeorico(a,b,pressureNormal,pInterferencia,a:0.1:b),'r',b:0.1:c,uTeorico(b,c,pInterferencia,0,b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','northwest','FontSize',16)
xlabel('r [mm]','FontSize',16)
ylabel('mm','FontSize',16)

%% Plots tension y malla deformada
% Deformed plot
figure
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% meshPlot(elements(1:nNodes1,:),nodes(1:nNodes1,:)+magnificationFactor*reshape(displacementsVector(1:nNodes1),[],nNodes1)','b','No');

% Von Mises Plot
Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));
figure; subplot(1,2,1); title('\sigma_{VM} [Mpa]','FontSize',14)
meshPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])','b','No');
bandPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])',vonMisesStress(1:nNodes1,:));
 subplot(1,2,2); title('\sigma_{VM} [Mpa]','FontSize',14)
meshPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])','b','No');
bandPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])',vonMisesStress(nElements1+1:end,:));

% toda la malla VM Plot
figure; title('Svm [MPa]')
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])',vonMisesStress);


% Sxx plot
Sxx = squeeze(elementStressExtrapolated(:,:,1)); 
figure; subplot(1,2,1); title('Int Sxx [Mpa]')
meshPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])','b','No');
bandPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])',Sxx(1:nNodes1,:));
 subplot(1,2,2); title('Out Sxx [Mpa]')
meshPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])','b','No');
bandPlot(elements2,nodes2+magnificationFactor*reshape(displacementsVector(convertNode2Dof(nNodes1+1:nNodes2+nNodes1,nDimensions)),nDimensions,[])',Sxx(nElements1+1:end,:));



