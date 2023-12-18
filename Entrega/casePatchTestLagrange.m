%% Case launcher
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

E = 200e3; nu = 0.3;
pressureNormal = 100;
a = 300; b = a + 108.1; c = b + 78.1; h = 600;
eEsf = b-a;
interferencia = 0.32; precond = 1e7;
nElementsZ = 20; nElementsR = 15; nElementsInRadius = 50; distorsion = 0;   
planeStrainFlag = 0; numbering = 'No'; % Yes/No
tol = 0.1; centroEsfera = [0 h];

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Axisymmetric';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              % Problem dimension

mallador

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

nodosSEsfSuperior = find(ismembertol(vecnorm(nodes-repmat(centroEsfera, size(nodes, 1), 1),2,2),a)); % nodos sobre la superficie interior de la semiesfera superior
nodosSEsfSuperior = nodosSEsfSuperior([2:end 1]);
pointLoadsArray = distributedLoad_Esfera(elementType,nodosSEsfSuperior,pointLoadsArray,nodes,pressureNormal,centroEsfera);

scatter(nodes([Lside nodosSEsfSuperior'],1),nodes([Lside nodosSEsfSuperior'],2),'filled','yellow')

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
pInterferenciaInt = 38.9;
pInterferenciaOut = 38.9;
pInterferencia = 38.9;


if planeStrainFlag
    E = E/(1-nu^2);
    nu = nu/(1-nu);
    fprintf("Caso plane strain \n")
else
    fprintf("Caso plane stress \n")
end

% pInterferencia=(E*interferencia/b) * (b^2-a^2)*(c^2-b^2) / (2*b^2*(c^2-a^2));
pInterferencia = abs(mean(mean(elementStressAtGaussPoints([90 301],:,1))))

C1= @(a,b,pInt,pOut) (a^2*pInt-b^2*pOut)/(b^2-a^2);
C2= @(a,b,pInt,pOut) (pInt-pOut)*a^2*b^2/(b^2-a^2);

sigmaR1=@(r)    C1(a,b,pressureNormal,pInterferencia)-C2(a,b,pressureNormal,pInterferencia)./r.^2;
sigmaTita1=@(r) C1(a,b,pressureNormal,pInterferencia)+C2(a,b,pressureNormal,pInterferencia)./r.^2;

sigmaR2=@(r)    C1(b,c,pInterferencia,0)-C2(b,c,pInterferencia,0)./r.^2;
sigmaTita2=@(r) C1(b,c,pInterferencia,0)+C2(b,c,pInterferencia,0)./r.^2;

nElementsCil = size(elements1,1);
figure; subplot(1,2,1); hold on; title('Srr'); grid
for iElements=1:nElementsCil
    scatter(nodes(elements(iElements,:),1),elementStressAtGaussPoints(iElements,:,1))
end
plot(a:0.1:b,sigmaR1(a:0.1:b),'r',b:0.1:c,sigmaR2(b:0.1:c),'b')

subplot(1,2,2); hold on; title('Stita'); grid
for iElements=1:nElementsCil
    scatter(nodes(elements(iElements,:),1),elementStressAtGaussPoints(iElements,:,2))
end
plot(a:0.1:b,sigmaTita1(a:0.1:b),'r',b:0.1:c,sigmaTita2(b:0.1:c),'b')

uTeorico= @(a,b,pInt,pOut,r) ((1-nu)/E).* C1(a,b,pInt,pOut) .* r + ((1+nu)/E).*  C2(a,b,pInt,pOut) ./r;

figure; hold on; title('u'); grid
for iElements=1:nElementsCil
    scatter(nodes(elements(iElements,:),1),displacementsMatrix(elements(iElements,:),1))
end
plot(a:0.1:b,uTeorico(a,b,pressureNormal,pInterferencia,a:0.1:b),'r',b:0.1:c,uTeorico(b,c,pInterferenciaOut,0,b:0.1:c),'b')

%% Solucion teorica semiesfera
% eEsf = 64.1;
% ri = a; ro = a+eEsf; pi = pressureNormal; po = 0;
% sigmaRR = @(r) po*ro.^3*(r.^3-ri.^3)./(r.^3*(ri.^3-ro.^3)) + pi*ri.^3*(ro.^3-r.^3)./(r.^3*(ri.^3-ro.^3));
% sigmaTT = @(r) po*ro.^3*(2*r.^3+ri.^3)./(2*r.^3*(ri.^3-ro.^3)) - pi*ri.^3*(2*r.^3+ro.^3)./(2*r.^3*(ri.^3-ro.^3));
% uEsfTeo = @(r) 3*(pi-po)*ri.^3*ro.^3./(4*E*r.^2*(ro.^3-ri.^3));
% sEsfVM = @(r) sqrt(sigmaRR(r).^2 + sigmaTT(r).^2 - sigmaRR(r).*sigmaTT(r));
% 
% figure('Name','S-Esf'); 
% subplot(1,4,1); hold on; title('Srr'); grid;
% plot(a:0.1:a+eEsf,sigmaRR(a:0.1:a+eEsf))
% subplot(1,4,2); hold on; title('Stt'); grid;
% plot(a:0.1:a+eEsf,sigmaTT(a:0.1:a+eEsf))
% subplot(1,4,3); hold on; title('Svm'); grid;
% plot(a:0.1:a+eEsf,sEsfVM(a:0.1:a+eEsf))
% subplot(1,4,4); hold on; title('u'); grid
% plot(a:0.1:a+eEsf,uEsfTeo(a:0.1:a+eEsf))

%% otra solucion, esta funciona PERO EN COMPRESION LRPM
% eEsf = 64.1;
% C1Esf = @(a,b,p) -p*b^3*(1-2*nu)/(E*(b^3-a^3));
% C2Esf = @(a,b,p) -p*a^3*b^3*(1+nu)/(2*E*(b^3-a^3));
% 
% simgaRR = @(r) ((1+nu)*C1Esf(a,a+eEsf,-pressureNormal)-2*(1-2*nu)*C2Esf(a,a+eEsf,-pressureNormal)./r.^3)*E/((1+nu)*1-2*nu);
% sigmaTT = @(r) ((1+nu)*C1Esf(a,a+eEsf,-pressureNormal)+(1-2*nu)*C2Esf(a,a+eEsf,-pressureNormal)./r.^3)*E/((1+nu)*1-2*nu);
% uEsfTeo = @(r) C1Esf(a,a+eEsf,-pressureNormal).*r + C2Esf(a,a+eEsf,-pressureNormal)./r.^2;
% sEsfVM = @(r) sqrt(sigmaRR(r).^2 + sigmaTT(r).^2 - sigmaRR(r).*sigmaTT(r));
% 
% figure('Name','S-Esf'); 
% subplot(1,4,1); hold on; title('Srr'); grid;
% plot(a:0.1:a+eEsf,sigmaRR(a:0.1:a+eEsf))
% subplot(1,4,2); hold on; title('Stt'); grid;
% plot(a:0.1:a+eEsf,sigmaTT(a:0.1:a+eEsf))
% subplot(1,4,3); hold on; title('Svm'); grid;
% plot(a:0.1:a+eEsf,sEsfVM(a:0.1:a+eEsf))
% subplot(1,4,4); hold on; title('u'); grid
% plot(a:0.1:a+eEsf,uEsfTeo(a:0.1:a+eEsf))

%% Plots tension y malla deformada
% Deformed plot
figure
meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% meshPlot(elements(1:nNodes1,:),nodes(1:nNodes1,:)+magnificationFactor*reshape(displacementsVector(1:nNodes1),[],nNodes1)','b','No');

% % Stresses plot
% figure
% meshPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)','b','No');
% bandPlot(elements,nodes+magnificationFactor*reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtGaussPoints(:,:,3)));
% title('Szz en los nodos [Mpa]')

% Von Mises Plot
Sxx = squeeze(elementStressExtrapolated(:,:,1))'; Syy = squeeze(elementStressExtrapolated(:,:,2))'; Szz = squeeze(elementStressExtrapolated(:,:,3))'; Szx = squeeze(elementStressExtrapolated(:,:,4))';
vonMisesStress = transpose(1/sqrt(2)*sqrt( (Sxx-Syy).^2 + (Syy-Szz).^2 + (Szz-Sxx).^2 + 6*Szx.^2 ));
figure; subplot(1,2,1); title('Int Svm [Mpa]')
meshPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])','b','No');
bandPlot(elements1,nodes1+magnificationFactor*reshape(displacementsVector(convertNode2Dof(1:nNodes1,nDimensions)),nDimensions,[])',vonMisesStress(1:nNodes1,:));
 subplot(1,2,2); title('Out Svm [Mpa]')
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



