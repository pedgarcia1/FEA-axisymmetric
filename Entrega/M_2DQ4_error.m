clear
clc
close all
format short g

load('caseFinalError.mat');
distanceTolerance = 0.01;

dofPerNode = 2;                    % Numero de grados de libertad por nodo
nodesPerElement = 4; % Numero de nodos por elemento
nElementalDof = dofPerNode*nodesPerElement;
numberOfElements = size(elements,1); 
nElements = numberOfElements; % Numero de elementos
numberOfNodes = size(nodes,1);  
nNodes = numberOfNodes; % Numero de nodos
totalNumberOfDof = dofPerNode*numberOfNodes;         % Numero de grados de libertad
nodeLabels = 1:numberOfNodes;

bc = boundaryConditionsArray;
R = pointLoadsArray;
C = constitutiveMatrix;   
stiffnessMatrix = K;
nodeDofs = reshape(1:totalNumberOfDof,dofPerNode,numberOfNodes)';

%% Gauss           
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
gaussPointPositions = [ -a  -a
                         a  -a
                         a   a
                        -a   a ];    
% Numero de puntos de Gauss
numberOfGaussPoints = size(gaussPointPositions,1);
gaussWeights = ones(numberOfGaussPoints,1);

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

% Solver
Dr = displacementsReducedVector;

% Reconstruccion
% D = zeros(totalNumberOfDof,1);
% D(isFree) = D(isFree) + Dr;
% displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector(1:end-nConstraints);
D = displacementsVector;

% Reacciones
% Rv = K(isFixed,isFree)*D(isFree);
% reacciones = nan(totalNumberOfDof,1);
% reacciones(isFixed) = Rv;
% reacciones = (reshape(reacciones,dofPerNode,[]))';

%% Recuperacion de tensiones en los puntos de Gauss
% stress = zeros(numberOfElements,nodesPerElement,3);
% 
% a = 1/sqrt(3);
% 
% stressPoints = [ -a  -a
%                  a  -a
%                  a   a
%                 -a   a ]; 
% 
% for iEle = 1:numberOfElements
%     elementNodes = nodes(elements(iEle,:),:);
%     for iNode = 1:nodesPerElement
%         % Punto de Gauss
%         ksi = stressPoints(iNode,1);
%         eta = stressPoints(iNode,2);  
%         % Derivadas de las funciones de forma respecto de ksi, eta
%         dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
%                   -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
%         % Derivadas de x,y, respecto de ksi, eta
%         jacobian = dN*elementNodes;                      
%         % Derivadas de las funciones de forma respecto de x,y.
%         dNxy = jacobian\dN;          % dNxy = inv(jacobian)*dN
% 
%         B = zeros(size(C,2),dofPerNode*nodesPerElement);
%         B(1,1:2:7) = dNxy(1,:);
%         B(2,2:2:8) = dNxy(2,:);
%         B(3,1:2:7) = dNxy(2,:);
%         B(3,2:2:8) = dNxy(1,:);
% 
%         eleDofs = nodeDofs(elements(iEle,:),:);
%         eleDofs = reshape(eleDofs',[],1);
%         stress(iEle,iNode,:) = C*B*D(eleDofs);
%     end
% end

stress = stressRecovery_2(elementType,elements,nodes,constitutiveMatrix,D);
% tension directa en los nodos

% [elementStressExtrapolated,elementStressAtGaussPoints] = stressRecext(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% stress = elementStressExtrapolated;


%% tensiones en el centro de cada elemento (superConvergencia)

stressCentro = zeros(numberOfElements,nodesPerElement,4);

for iEle = 1:numberOfElements
    elementNodes = nodes(elements(iEle,:),:);
    for iNode = 1:nodesPerElement
        % Punto de Gauss
        ksi = 0;
        eta = 0;  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jacobian = dN*elementNodes;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jacobian\dN;          % dNxy = inv(jacobian)*dN
        
        % B = zeros(size(C,2),dofPerNode*nodesPerElement);
        % B(1,1:2:7) = dNxy(1,:);
        % B(2,2:2:8) = dNxy(2,:);
        % B(3,1:2:7) = dNxy(2,:);
        % B(3,2:2:8) = dNxy(1,:);

        N4 = 0.25*(1 - ksi)*(1 + eta);
        N3 = 0.25*(1 + ksi)*(1 + eta);
        N2 = 0.25*(1 + ksi)*(1 - eta);
        N1 = 0.25*(1 - ksi)*(1 - eta);

        Ni = [N1 N2 N3 N4];

        r = Ni*elementNodes(:,1);

        % jacoVec(elements(iele,iNode)) = det(jacobian);

        B(1,1:2:nElementalDof-1) = dNxy(1,:);
        B(2,1:2:nElementalDof-1) =  Ni/r; % N(xi,eta)_i/r(xi,eta) para cada pto de gauss
        B(3,2:2:nElementalDof) = dNxy(2,:);
        B(4,1:2:nElementalDof-1) = dNxy(2,:);
        B(4,2:2:nElementalDof) = dNxy(1,:);

        eleDofs = nodeDofs(elements(iEle,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressCentro(iEle,iNode,:) = C*B*D(eleDofs);
    end
end

%% promediado de tensiones en los nodos

avgStress = zeros(numberOfNodes,4);

for iNode = 1:numberOfNodes
    [I,J] = find(elements == iNode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(iNode,:) = avgStress(iNode,:) + squeeze(stressCentro(I(ishare),J(ishare),:))';
%         si no todos los elementos tienen la misma area, deberian tener
%         pesos distintos!! eso me puede distorsionar los valores.
    end
    avgStress(iNode,:) = avgStress(iNode,:) / nShare;
end

%% calculo de eta_el, e2, U2
invC = C\eye(4);
eta_el = zeros(numberOfElements,1);
e2_el = zeros(numberOfElements,1);
U2_el = zeros(numberOfElements,1);

for iEle = 1:numberOfElements
    elementNodes = nodes(elements(iEle,:),:);
    for ipg = 1:numberOfGaussPoints
        % Punto de Gauss
        ksi = gaussPointPositions(ipg,1);
        eta = gaussPointPositions(ipg,2);
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jacobian = dN*elementNodes;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jacobian\dN;          % dNxy = inv(jacobian)*dN
        
        % funciones de forma
        N4 = 0.25*(1 - ksi)*(1 + eta);
        N3 = 0.25*(1 + ksi)*(1 + eta);
        N2 = 0.25*(1 + ksi)*(1 - eta);
        N1 = 0.25*(1 - ksi)*(1 - eta);
        N = [ N1 N2 N3 N4 ];
        eleStress =  squeeze(stress(iEle,ipg,:));            % tensiones "directas"::: em los nodos
        starStress = ( N * avgStress(elements(iEle,:),:) )'; % tensiones mejoradas::: en el centro tenes superconvergente del Q4.. 
%                                 promedias los centro de cada elemento a cada nodo y con eso tenes las tensiones mejoradas
        
        e2_el(iEle) = e2_el(iEle) + (starStress - eleStress)' * ... 
                    invC * (starStress - eleStress) * gaussWeights(ipg) * det(jacobian);
                
        U2_el(iEle) = U2_el(iEle) + eleStress' * invC * eleStress * ...
                      gaussWeights(ipg) * det(jacobian);
    end
    
    eta_el(iEle) = sqrt( e2_el(iEle) / (e2_el(iEle) + U2_el(iEle)) );
    
end

etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) );


%% Configuracion deformada
D = (reshape(D,dofPerNode,[]))';
nodePosition = nodes + D(:,1:2);

%Graficacion
limites = [-100 0];

labelNames = cellstr(num2str((1:size(nodes,1))'));

figure(1)
bandploterror(elements,nodePosition,stressCentro(:,:,1),limites,'k');
title('Tensiones en el centro de cada elemento.')
axis equal tight

figure(2)
scalarbandplot(elements,nodePosition,avgStress(:,1),limites,'k',[],'interp');
title('Tensiones promediadas.')
axis equal tight

figure(3)
scalarbandplot(elements,nodePosition,eta_el,[],'k',[],'flat');
title(' $\eta$ (Error  ZZ)','Interpreter','latex')
axis equal tight

figure(4)
meshPlot(elements,nodePosition,'k');
% text(nodePosition(:,1),nodePosition(:,2), labelNames)
% text(nodePosition(:,1),nodePosition(:,2), [])
title(' $\eta$ (Error  ZZ)','Interpreter','latex')
axis equal tight

figure(5)
bandploterror(elements,nodes,stress(:,:,1),limites,'k')
title('Tension directa en los nodos')
axis equal tight


