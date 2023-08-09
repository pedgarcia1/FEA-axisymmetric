clc; clear; close all; 

keyPlots=0;
keyPlots2=0;

radio=100;
LargoXIn=70;
LargoXOut=20;
interferencia=0.1;


nElemEnX=20;
nElemEnXOut=10;
nElemEnY=10;
AlturaY=100;


p0=200*radio*AlturaY/nElemEnY;

nDimensions=2; 
eleType='Q4';
precondCT=1e7;
E=200000;
poisson=0.3;
%% Mesh generation
 [meshInfo.elements,meshInfo.nodes,meshInfo.vertexNodes,meshInfo.sideNodes]=quadrilateralDomainMeshGenerator(eleType,'Straight',LargoXIn,AlturaY,nElemEnX,nElemEnY,0,0);  
 meshInfo.nodes(:,1)=meshInfo.nodes(:,1)+radio; %Axisimetrico respecto de X

 [meshInfoOut.elements,meshInfoOut.nodes,meshInfoOut.vertexNodes,meshInfoOut.sideNodes]=quadrilateralDomainMeshGenerator(eleType,'Straight',LargoXOut,AlturaY,nElemEnXOut,nElemEnY,0,0);  
 meshInfoOut.nodes(:,1)=meshInfoOut.nodes(:,1)+radio+LargoXIn; %Los hacemos coincidir aunque sin interferencia



nElements1=size(meshInfo.elements,1);    %Number of elements
nElements2=size(meshInfoOut.elements,1);
nElemTotal=nElements1+nElements2;
nNodes1=size(meshInfo.nodes,1);          %Number of nodes
nNodes2=size(meshInfoOut.nodes,1);          %Number of nodes
nNodesTotal=nNodes1+nNodes2;
nTotalDof=nNodesTotal*nDimensions;           %Number of total dofs

if keyPlots
    figure; hold on
    nodosOutPlot=meshInfoOut.nodes;
    nodosOutPlot(:,1)=nodosOutPlot(:,1)+radio;
    meshPlot(meshInfo.elements,meshInfo.nodes,'b','Yes');
    meshPlot(meshInfoOut.elements,nodosOutPlot,'b','Yes');
    axis([0 max(nodosOutPlot(:,1)) 0 max(meshInfoOut.nodes(:,2))])
    
    figure; hold on; title('Malla Inicial')
    nodosOutPlot=meshInfoOut.nodes;
    nodosOutPlot(:,1)=nodosOutPlot(:,1)-interferencia;
    meshPlot(meshInfo.elements,meshInfo.nodes,'b','Yes');
    meshPlot(meshInfoOut.elements,nodosOutPlot,'b','Yes');
    axis([0 max(nodosOutPlot(:,1)) 0 max(meshInfoOut.nodes(:,2))])
end

nDofTotal=2*nNodesTotal;
nodeDofs = reshape(1:nDofTotal, 2, [])'; 

%% Buscamos los nodos que se van a compartir
tolerancia = 1e-6;

constraintsRelations = [];

for i = 1:size(meshInfo.nodes, 1)
    nodo_malla1 = meshInfo.nodes(i, :);    
    for j = 1:size(meshInfoOut.nodes, 1)
        nodo_malla2 = meshInfoOut.nodes(j, :);
        
        % Calcula la distancia entre los nodos en cada coordenada
        distancia = norm(nodo_malla1 - nodo_malla2);
        
        % Compara la distancia con la tolerancia
        if distancia < tolerancia
            constraintsRelations = [constraintsRelations; i, j];
        end
    end
end

constraintsRelations(:,2)=constraintsRelations(:,2)+nNodes1; %CHEQUEAR ESTO

meshInfoOut.nodes(:,1)=meshInfoOut.nodes(:,1)-interferencia; 
% Segun lo
% que probe es lo mismo quien se lleve la deformacion incial, 
% duda, hay que ponerla o el mismo sistema la impone?
% meshInfoOut.nodes(:,2)=meshInfoOut.nodes(:,2)+interferencia; % 

%% Armamos K

[constitutive] = constitutiveIsotropicMatrix('Axisymmetric',E,poisson);

K1=assembleStiffnessMatrixAxisimetric(meshInfo.nodes,meshInfo.elements,constitutive,eleType,'x');
K2=assembleStiffnessMatrixAxisimetric(meshInfoOut.nodes,meshInfoOut.elements,constitutive,eleType,'x');

%% Aplicamos Cargas y BC

% Boundary conditions
boundaryConditionsArrayIn = false(nNodes1,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArrayIn(meshInfo.sideNodes([1 3],:),2) = true;

boundaryConditionsArrayOut = false(nNodes2,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArrayOut(meshInfoOut.sideNodes([1 3],:),2) = true;



%% Aplicamos cargas
pointLoadsArrayIN = zeros(nNodes1,nDimensions);            % Point load nodal value for each direction

pointLoadsArrayIN(meshInfo.sideNodes(4,meshInfo.sideNodes(4,:)>0),1)=p0;
pointLoadsArrayIN(meshInfo.vertexNodes([1 3]),1)=p0/2;


pointLoadsArrayOut = zeros(nNodes2,nDimensions);  

% SideNodes tiene abajo derecha arriba izquierda

%% PatchTest

% pointLoadsArrayIN(meshInfo.sideNodes(2,meshInfo.sideNodes(2,:)>0),1)=-p0*(radio+LargoXIn)/radio;
% pointLoadsArrayIN(meshInfo.vertexNodes([2 4]),1)=-p0/2*(radio+LargoXIn)/radio;

%% Hacemos matriz de CT


[nCT,CTMatrix] = getCTzunchado(nodeDofs,constraintsRelations,precondCT,nDofTotal);

%% Solver

% Matrix reduction
isFixedFull = reshape([boundaryConditionsArrayIn; boundaryConditionsArrayOut]',1,[])';
isFreeFull = ~isFixedFull;


assert(K1(1,1)>0.1,'Diagonal Neg')
assert(K2(1,1)>0.1,'Diagonal Neg')

% Loads Vector rearrangement
loadsVectorFull = reshape([pointLoadsArrayIN; pointLoadsArrayOut]',1,[])';

%% Ensamblamos 

K=[K1 zeros(size(K1,2),size(K2,1));zeros(size(K2,1),size(K1,2)) K2];

KFull=[K CTMatrix';CTMatrix zeros(nCT)];

R=[loadsVectorFull(isFreeFull);ones(nCT,1).*interferencia.*precondCT];

isFreeFullCT=[isFreeFull; true(nCT,1)];

% Equation solving
displacementsReducedVector = KFull(isFreeFullCT,isFreeFullCT)\R;

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFreeFull) = displacementsVector(isFreeFull) + displacementsReducedVector(1:end-nCT);

if keyPlots 
    figure;hold on;title('Malla Deformada')
    meshPlot(meshInfo.elements,meshInfo.nodes+reshape(displacementsVector(1:nNodes1*2),nDimensions,nNodes1)','b','Yes');
    meshPlot(meshInfoOut.elements,meshInfoOut.nodes+reshape(displacementsVector(nNodes1*2+1:end),nDimensions,nNodes2)','b','Yes');
    if keyPlots2
        for iNode=1:nNodes

            if boundaryConditionsArray(iNode,1) && boundaryConditionsArray(iNode,2)
                scatter(meshInfo.nodes(iNode,1),meshInfo.nodes(iNode,2),'black', 'filled')
            elseif boundaryConditionsArray(iNode,1)
                scatter(meshInfo.nodes(iNode,1),meshInfo.nodes(iNode,2),'blue', 'filled')
            elseif boundaryConditionsArray(iNode,2)
                scatter(meshInfo.nodes(iNode,1),meshInfo.nodes(iNode,2),'yellow', 'filled')
            end
            legend('Restringido en: Black ambos, blue X, yellow Y')

            quiver(meshInfo.nodes(iNode,1),meshInfo.nodes(iNode,2), pointLoadsArrayIN(iNode, 1)./10e6, pointLoadsArrayIN(iNode, 2)./10e6, 'color', 'k', 'linewidth', 1.5, 'maxheadsize', 0.5)
        end
    end
end

%% Calculamos Tensiones

[elementStressAtGaussPoints]=stressRecovery(eleType,meshInfo.elements,meshInfo.nodes,constitutive,displacementsVector(1:nNodes1*2));
[elementStressAtGaussPointsOut]=stressRecovery(eleType,meshInfoOut.elements,meshInfoOut.nodes,constitutive,displacementsVector(nNodes1*2+1:end));
    
% figure;hold on;title('Tensiones r')
% plotColo(meshInfo.nodes,meshInfo.elements,elementStressAtGaussPoints(:,:,1)); %% Estas Tensiones dan perfecto (Tensiones en r)
% plotColo(meshInfoOut.nodes,meshInfoOut.elements,elementStressAtGaussPointsOut(:,:,1)) %% Estas Tensiones dan perfecto (Tensiones en r)
% 
% figure;hold on;title('Tensiones tita')
% plotColo(meshInfo.nodes,meshInfo.elements,elementStressAtGaussPoints(:,:,2)); %% Estas Tensiones dan perfecto (Tensiones en r)
% plotColo(meshInfoOut.nodes,meshInfoOut.elements,elementStressAtGaussPointsOut(:,:,2)) %% Estas Tensiones dan perfecto (Tensiones en r)

figure;
subplot(1,2,1);hold on;title('Tensiones VM IN')
plotColo(meshInfo.nodes,meshInfo.elements,elementStressAtGaussPoints(:,:,5)); %% Estas Tensiones dan perfecto (Tensiones en r)
subplot(1,2,2);hold on;title('Tensiones VM Out')
plotColo(meshInfoOut.nodes,meshInfoOut.elements,elementStressAtGaussPointsOut(:,:,5)) %% Estas Tensiones dan perfecto (Tensiones en r)


%% Plot desplazamientos

% despla=reshape(displacementsVector,2,[])';
% figure; hold on

% for iNodes=1:nNodesTotal
%     scatter(meshInfo.nodes(iNodes,1),despla(iNodes,1))
% end

%% comparamos 
%  CheckResults
% 
% fprintf('Solucion teorica %f %f \n',[ur(a) ur(b)])
% despla=reshape(displacementsVector,2,[])';
% solFea=max(abs(despla));
% fprintf('Solucion Fea     %.6f %.6f \n',[solFea(1) solFea(2)])

%% Solucion teorica
a=radio;b=radio+LargoXIn;c=b+LargoXOut;

p=(E*interferencia/b) * (b^2-a^2)*(c^2-b^2) / (2*b^2*(c^2-a^2));

C1= @(a,b,pInt,pOut) (a^2*pInt-b^2*pOut)/(b^2-a^2);
C2= @(a,b,pInt,pOut) (pInt-pOut)*a^2*b^2/(b^2-a^2);

sigmaRIn=@(r)    C1(a,b,200,p)-C2(a,b,200,p)./r.^2;
sigmaTitaIn=@(r) C1(a,b,200,p)+C2(a,b,200,p)./r.^2;

sigmaROut=@(r)    C1(b,c,p,0)-C2(b,c,p,0)./r.^2;
sigmaTitaOut=@(r) C1(b,c,p,0)+C2(b,c,p,0)./r.^2;

%% Ploteamos tensiones FEA
% Sigma R IN
subplot(2,2,1);hold on; title('Comparativa Sigma R IN')

for iElements=1:nElements1
    scatter(meshInfo.nodes(meshInfo.elements(iElements,:),1),elementStressAtGaussPoints(iElements,:,1))
end
plot(radio:0.1:LargoXIn+radio,sigmaRIn(radio:0.1:LargoXIn+radio),'r');

%Sigma Tita IN
subplot(2,2,2);hold on; title('Comparativa Sigma Tita IN')

for iElements=1:nElements1
    scatter(meshInfo.nodes(meshInfo.elements(iElements,:),1),elementStressAtGaussPoints(iElements,:,2))
end
plot(radio:0.1:LargoXIn+radio,sigmaTitaIn(radio:0.1:LargoXIn+radio),'r');


% Sigma R Out
subplot(2,2,3);hold on; title('Comparativa Sigma R Out')

for iElements=1:nElements2
    scatter(meshInfoOut.nodes(meshInfoOut.elements(iElements,:),1),elementStressAtGaussPointsOut(iElements,:,1))
end
plot(LargoXIn+radio:0.1:LargoXOut+LargoXIn+radio,sigmaROut(LargoXIn+radio:0.1:LargoXOut+LargoXIn+radio),'r');

%Sigma Tita Out
subplot(2,2,4);hold on; title('Comparativa Sigma Tita Out')

for iElements=1:nElements2
    scatter(meshInfoOut.nodes(meshInfoOut.elements(iElements,:),1),elementStressAtGaussPointsOut(iElements,:,2))
end
plot(LargoXIn+radio:0.1:LargoXOut+LargoXIn+radio,sigmaTitaOut(LargoXIn+radio:0.1:LargoXOut+LargoXIn+radio),'r');
