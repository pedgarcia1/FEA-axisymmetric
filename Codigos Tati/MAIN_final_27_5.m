clc; clear; close all; 

keyPlots=0;
keyPlots2=0;
radio=100;
LargoX=400;
nElemEnY=20;
AlturaY=200;


p0=200*radio*AlturaY/nElemEnY;

nDimensions=2; 
eleType='Q4';


%% Mesh generation
 [meshInfo.elements,meshInfo.nodes,meshInfo.vertexNodes,meshInfo.sideNodes]=quadrilateralDomainMeshGenerator(eleType,'Straight',LargoX,AlturaY,40,nElemEnY,0,5);
%  
 meshInfo.nodes(:,1)=meshInfo.nodes(:,1)+radio; %Axisimetrico respecto de X

nElements=size(meshInfo.elements,1);    %Number of elements
nNodes=size(meshInfo.nodes,1);          %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

if keyPlots
    meshPlot(meshInfo.elements,meshInfo.nodes,'b','Yes');
    axis([0 max(meshInfo.nodes(:,1)) 0 max(meshInfo.nodes(:,2))])
end
%% Armamos K

[constitutive] = constitutiveIsotropicMatrix('Axisymmetric',200000,0.3);

K=assembleStiffnessMatrixAxisimetric(meshInfo.nodes,meshInfo.elements,constitutive,eleType,'x');

%% Aplicamos Cargas y BC

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(meshInfo.sideNodes([1 3],:),2) = true;




%% Primero tenemos que identificar los nodos de los bordes ya que tienen la mitad de la carga

% numerosComunes_1_4 = intersect(meshInfo.sideNodes(1,meshInfo.sideNodes(1,:)>0), meshInfo.sideNodes(4,meshInfo.sideNodes(4,:)>0));
% 
% numerosComunes_4_3 = intersect(meshInfo.sideNodes(4,meshInfo.sideNodes(4,:)>0), meshInfo.sideNodes(3,meshInfo.sideNodes(3,:)>0));

%% aplicamos cargas
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
pointLoadsArray(meshInfo.sideNodes(4,meshInfo.sideNodes(4,:)>0),1)=p0;
pointLoadsArray(meshInfo.sideNodes(2,meshInfo.sideNodes(2,:)>0),1)=-p0*(radio+LargoX)/radio;

% SideNodes tiene abajo derecha arriba izquierda

% pointLoadsArray([numerosComunes_1_4 numerosComunes_4_3],1)=p0/2;

pointLoadsArray(meshInfo.vertexNodes([1 3]),1)=p0/2;
pointLoadsArray(meshInfo.vertexNodes([2 4]),1)=-p0/2*(radio+LargoX)/radio;
%% Solver

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% assert(min(sort(abs(eig(K(isFree,isFree)))))>0.1,'Aval 0')
assert(K(1,1)>0.1,'Diagonal Neg')

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = K(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;

if keyPlots && keyPlots2
    figure
    meshPlot(meshInfo.elements,meshInfo.nodes+reshape(displacementsVector,nDimensions,nNodes)','b','Yes');
    hold on
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

            quiver(meshInfo.nodes(iNode,1),meshInfo.nodes(iNode,2), pointLoadsArray(iNode, 1)./10e6, pointLoadsArray(iNode, 2)./10e6, 'color', 'k', 'linewidth', 1.5, 'maxheadsize', 0.5)
        end
    end
end

%% Calculamos Tensiones

[elementStressAtGaussPoints]=stressRecovery(eleType,meshInfo.elements,meshInfo.nodes,constitutive,displacementsVector);

maximos=[max(elementStressAtGaussPoints(:,:,1)) min(elementStressAtGaussPoints(:,:,1)); max(elementStressAtGaussPoints(:,:,2)) min(elementStressAtGaussPoints(:,:,2)); max(elementStressAtGaussPoints(:,:,3)) min(elementStressAtGaussPoints(:,:,3)); max(elementStressAtGaussPoints(:,:,4)) min(elementStressAtGaussPoints(:,:,4))]; 

figure;plotColo(meshInfo.nodes,meshInfo.elements,elementStressAtGaussPoints(:,:,1));title('Tensiones r') %% Estas Tensiones dan perfecto (Tensiones en r)

figure;plotColo(meshInfo.nodes,meshInfo.elements,elementStressAtGaussPoints(:,:,2));title('Tensiones tita') %% Estas Tensiones no dan (tensiones en tita)

%% Plot desplazamientos

despla=reshape(displacementsVector,2,[])';
figure; hold on

for iNodes=1:nNodes
    scatter(meshInfo.nodes(iNodes,1),despla(iNodes,1))
end

%% comparamos 
 CheckResults

fprintf('Solucion teorica %f %f \n',[ur(a) ur(b)])
despla=reshape(displacementsVector,2,[])';
solFea=max(abs(despla));
fprintf('Solucion Fea     %.6f %.6f \n',[solFea(1) solFea(2)])
