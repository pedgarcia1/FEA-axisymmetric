% convergencia de malla
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

addpath('Convergencia de Malla\')
archivos = dir("Convergencia de Malla\");
for i = 3:length(archivos)
    fileNames{i-2} = archivos(i).name;
end

nodConcPos = [300.0000 -500.0000;
     275.2124  119.4073;
    297.5770  -1105.2388;
    0.1      -1434.3929];

nCasos = round(length(fileNames)/3);
for iCaso = 2:nCasos+1

mshTxt =  sprintf('final_ref_%d.dat',iCaso);
loadsTxt = sprintf('finalLoads_ref_%d.csv',iCaso);
constraintsTxt = sprintf('finalConstraints_ref_%d.csv',iCaso);

%%    
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
msh = impq4(mshTxt);
msh = impq4('final_ref_7.dat')
nodes = msh.nodes; elements = msh.elements;
size(elements,1)
loads = impLoadsNX(loadsTxt);
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

msh.c = impConstraintsNX(constraintsTxt);

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
% [elementStressAtNodes]=stressRecovery_2(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

%% RESULT PLOTS
magnificationFactor=1;
displacementsMatrix = reshape(displacementsVector,nDimensions,nNodes)';  %mm

%% Von Mises Plot
Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));

%% Process stress and convergencia plot

for iConcNode= 1:size(nodConcPos,1)

    distances = sqrt((nodes(:,1) - nodConcPos(iConcNode,1)).^2 + (nodes(:,2) - nodConcPos(iConcNode,2)).^2 ) ;
    [~, closestIndex] = min(distances); stressInterpolationFactor(1) = 0.99;
    closestNode = nodes(closestIndex,:);

    if iConcNode == 1
        resultados(iCaso,iConcNode) = norm(displacementsMatrix(closestIndex,:));
    else
    log = closestIndex;
    for i = 1:size(log,1)
        eleLog = ismember(elements,log(i));
        iEles = find(any(eleLog')');
        for n = 1:size(iEles)
            auxStress(n) = vonMisesStress(iEles(n),eleLog(iEles(n),:)');
            auxStress(n) = auxStress(n)*stressInterpolationFactor(1);
        end
        resultados(iCaso,iConcNode) = mean(auxStress);
        auxStress = [];
    end
    end

end

resultados(iCaso,size(nodConcPos,1) +1) = nElements;

end
%%
mallasElementos = resultados(:,5)';

fConvergencia = figure('Name','Convergencia de Malla');
hold on
xlabel('Numero de Elementos')
grid
% yyaxis left
ylabel('\sigma_{VM} [MPa]')
plot(mallasElementos,resultados(:,2)')
plot(mallasElementos,resultados(:,3)')
plot(mallasElementos,resultados(:,4)')
% yyaxis right
% ylabel('u [mm]')
% plot(mallasElementos,resultados(:,1)')

legend('C1','C2','C3')

%% calculo error

e = [];
for i = 1:size(resultados,1)-1

    phi(1) = resultados(i,3);
    phi(2) = resultados(i+1,3);
    h(1) = (resultados(i,5)).^(-2/2);
    h(2) = (resultados(i+1,5)).^(-2/2);
    e(i) = abs(calcErr(phi,h));

end

figure; plot(resultados(2:end,5),e)