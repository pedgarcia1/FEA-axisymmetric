%% Case launcher
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 108.1; c = b + 78.1; h = 1200;
eEsf = b-a;
interferencia = 0.31; precond = 1e7;
nElementsZ = 20; nElementsR = 15; nElementsInRadius = 50; distorsion = 0;   
planeStrainFlag = 0; numbering = 'No'; % Yes/No
tol = 0.1; centroEsfera = [0 h];

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

% import NX mesh
msh = impq4('final.dat');
nodes = msh.nodes; elements = msh.elements;
loads = impLoadsNX('finalLoads.csv');
%

nElements=size(elements,1);    %Number of elements
nNodes=size(nodes,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(msh.bcond.fijox,1) = true;
boundaryConditionsArray(msh.bcond.fijoz,2) = true;

figure; hold on; meshPlot(elements,nodes,'b',numbering); scatH = [];
scatH(1) = scatter(nodes(boundaryConditionsArray(:,1),1),nodes(boundaryConditionsArray(:,1),2),'filled','red');
scatH(end+1) = scatter(nodes(boundaryConditionsArray(:,2),1),nodes(boundaryConditionsArray(:,2),2),'filled','green');

% Load definition
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
pointLoadsArray(loads.NodeID,:) = [loads.R,loads.Z];

scatH(end+1) = scatter(nodes(pointLoadsArray(:,1) ~= 0 | pointLoadsArray(:,2) ~= 0,1),nodes(pointLoadsArray(:,1) ~= 0 | pointLoadsArray(:,2) ~= 0,2),'filled','yellow');

%% Solver

% Constraints
% - u1 + u2 = interferencia
% - sideNodes1(2,:) + sideNodes2(4,:) = intereferencia
% nConstraints = length(sideNodes(2,:));

msh.c = impConstraintsNX('finalConstraints.csv');

maxR = max(msh.c.RCoord);
intCilNodes = msh.c.NodeID(msh.c.RCoord == maxR);
% intCilNodes = intCilNodes(2:end);
minR = min(msh.c.RCoord);
outCilNodes = msh.c.NodeID(msh.c.RCoord == minR);

nConstraints = length(intCilNodes); 
C = zeros(nConstraints,nTotalDof);
% rightSideDOF = convertNode2Dof(sideNodes(2,:),nDimensions);
% leftSideDOF = convertNode2Dof(sideNodes2(4,:),nDimensions);
rightSideDOF = convertNode2Dof(intCilNodes',nDimensions);
leftSideDOF = convertNode2Dof(outCilNodes',nDimensions);
n = 1;
for i = 1:2:length(rightSideDOF)
    C(n,rightSideDOF(i)) = -1*precond; % u1
    nod1 = (rightSideDOF(i) + 1)/nDimensions;
    [~,I] = sort(sqrt( (nodes(nod1,1) - nodes(outCilNodes,1)).^2 + (nodes(nod1,2) - nodes(outCilNodes,2)).^2 ),'ascend');
    nod2 = outCilNodes(I(1));
    nod2dof = convertNode2Dof(nod2,nDimensions);
    C(n,nod2dof(1)) = 1*precond; % u2
    n = n + 1;
end

scatH(end+1) = scatter(nodes(intCilNodes,1),nodes(intCilNodes,2),'black','filled');
scatH(end+1) = scatter(nodes(outCilNodes,1),nodes(outCilNodes,2),'blue','filled');
legend(scatH,{'fijo R','fijo Z','Pint','C int','C out'},'Location','bestoutside')

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

% pInterferencia = abs(mean(mean(elementStressAtGaussPoints([90 301],:,1))))
outCilElements = find(any(ismember(elements,outCilNodes)'));
intCilElements = find(any(ismember(elements,intCilNodes)'));
pInterferenciaInt = abs(mean(mean(elementStressAtGaussPoints(intCilElements,:,1))));
pInterferenciaOut = abs(mean(mean(elementStressAtGaussPoints(outCilElements,:,1))));
pInterferencia = pInterferenciaInt;

C1= @(a,b,pInt,pOut) (a^2*pInt-b^2*pOut)/(b^2-a^2);
C2= @(a,b,pInt,pOut) (pInt-pOut)*a^2*b^2/(b^2-a^2);

sigmaR1=@(r)    C1(a,b,pressureNormal,pInterferencia)-C2(a,b,pressureNormal,pInterferencia)./r.^2;
sigmaTita1=@(r) C1(a,b,pressureNormal,pInterferencia)+C2(a,b,pressureNormal,pInterferencia)./r.^2;
sigma_z1 = pressureNormal*a^2./(b.^2-a^2);

sigmaR2=@(r)    C1(b,c,pInterferencia,0)-C2(b,c,pInterferencia,0)./r.^2;
sigmaTita2=@(r) C1(b,c,pInterferencia,0)+C2(b,c,pInterferencia,0)./r.^2;

sigmaVM1=@(r) sqrt(0.5*((sigmaR1(r)-sigmaTita1(r)).^2 + (sigmaTita1(r)-sigma_z1).^2 + (sigma_z1-sigmaR1(r)).^2 ));
sigmaVM2=@(r) sqrt(0.5*((sigmaR2(r)-sigmaTita2(r)).^2 + (sigmaTita2(r)-0).^2 + (0-sigmaR2(r)).^2 ));

resultadosNX = importResultadosNX('finalResultadosNX.csv');

figure; subplot(1,2,1); hold on; title('\sigma_{r}','Interpreter','tex'); grid
% zPlot = -500;
zPlot = -500;
log = find(abs(nodes(:,2) - zPlot) < 10);
for i = 1:size(log,1)
    eleLog = ismember(elements,log(i));
    iEles = find(any(eleLog')');
    for n = 1:size(iEles)
        auxStress(n) = elementStressExtrapolated(iEles(n),eleLog(iEles(n),:)',1);
        try
            auxNXStress(n) = resultadosNX.RR(resultadosNX.ElemID == iEles(n) & resultadosNX.NodeID == log(i));
        catch
            auxNXStress = [];
        end
    end
    h(1) = scatter(nodes(log(i),1),mean(auxStress));
    h(2) = scatter(nodes(log(i),1),mean(auxNXStress),'x');
    auxStress = [];
end

aux = plot(a:0.1:b,sigmaR1(a:0.1:b),'r',b:0.1:c,sigmaR2(b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','northwest','FontSize',16)
xlabel('r [mm]','FontSize',16); ylabel('MPa','FontSize',16);
ac = gca();
ac.XAxis.FontSize = 16;
ac.YAxis.FontSize = 16;

subplot(1,2,2); hold on; title('\sigma_{\theta}','Interpreter','tex'); grid
log = find(abs(nodes(:,2) - zPlot) < 10);
for i = 1:size(log,1)
    eleLog = ismember(elements,log(i));
    iEles = find(any(eleLog')');
    for n = 1:size(iEles)
        auxStress(n) = elementStressExtrapolated(iEles(n),eleLog(iEles(n),:)',2);
        try
            auxNXStress(n) = resultadosNX.TT(resultadosNX.ElemID == iEles(n) & resultadosNX.NodeID == log(i));
        catch
            auxNXStress = [];
        end
    end
    h(1) = scatter(nodes(log(i),1),mean(auxStress));
    h(2) = scatter(nodes(log(i),1),mean(auxNXStress),'x');
    auxStress = [];
end

aux = plot(a:0.1:b,sigmaTita1(a:0.1:b),'r',b:0.1:c,sigmaTita2(b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','northwest','FontSize',16)
xlabel('r [mm]','FontSize',16); ylabel('MPa','FontSize',16);
ac = gca();
ac.XAxis.FontSize = 16;
ac.YAxis.FontSize = 16;

uTeorico= @(a,b,pInt,pOut,r) ((1-nu)/E).* C1(a,b,pInt,pOut) .* r + ((1+nu)/E).*  C2(a,b,pInt,pOut) ./r;
%%
figure; hold on; title('u'); grid
log = find(abs(nodes(:,2) - zPlot) < 10);
aux = scatter(nodes(log,1),displacementsMatrix(log,1));
h(1) = aux(1);
aux = plot(a:0.1:b,uTeorico(a,b,pressureNormal,pInterferencia,a:0.1:b),'r',b:0.1:c,uTeorico(b,c,pInterferencia,0,b:0.1:c),'b');
h(2) = aux(1);
legend(h,{'FEA','Sol. Teorica'},'Location','northwest','FontSize',16)

%% Deformed plot
figure; 
% title(sprintf('Deformed Plot MF: %d',magnificationFactor));
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% meshPlot(elements(1:nNodes1,:),nodes(1:nNodes1,:)+magnificationFactor*reshape(displacementsVector(1:nNodes1),[],nNodes1)','b','No');
% xlim([0 700])

%% Von Mises Plot
Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));

% toda la malla VM Plot
figure; title('\sigma_{VM} [MPa]','Interpreter','tex','FontSize',16)
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])','b','No');
bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,[])',vonMisesStress);
fallaNodes = vonMisesStress > 250*ones(nElements,4);
fallaNodes = elements(fallaNodes);
scatter(nodes(fallaNodes,1),nodes(fallaNodes,2),'red','filled');

%% Teorico VM Plot 

figure; hold on; title('\sigma_{VM}','Interpreter','tex','FontSize',16); grid
log = find(abs(nodes(:,2) - zPlot) < 10);
for i = 1:size(log,1)
    eleLog = ismember(elements,log(i));
    iEles = find(any(eleLog')');
    for n = 1:size(iEles)
        auxStress(n) = vonMisesStress(iEles(n),eleLog(iEles(n),:)');
        try
            auxNXStress(n) = resultadosNX.VonMises(resultadosNX.ElemID == iEles(n) & resultadosNX.NodeID == log(i));
        catch
            auxNXStress = [];
        end
    end
    h(1) = scatter(nodes(log(i),1),mean(auxStress));
    h(2) = scatter(nodes(log(i),1),mean(auxNXStress),'x');
    auxStress = [];
end

xlabel('r [mm]','FontSize',16); ylabel('MPa','FontSize',16);
aux = plot(a:0.1:b,sigmaVM1(a:0.1:b),'r',b:0.1:c,sigmaVM2(b:0.1:c),'b');
h(3) = aux(1);
legend(h,{'FEA','NX','Sol. Teorica'},'Location','northwest','FontSize',16)
ac = gca();
ac.XAxis.FontSize = 16;
ac.YAxis.FontSize = 16;


