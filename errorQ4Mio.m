%% mi error Q4
load('testCaseFinalError.mat');
C = constitutiveMatrix;
%% calculo de eta_el, e2, U2
invC = C\eye(4);
eta_el = zeros(nElements,1);
e2_el = zeros(nElements,1);
U2_el = zeros(nElements,1);

[gaussPointsLocation,gaussPointsWeight] = getGaussPoints('Quadrilateral',4);
npg=size(gaussPointsWeight,1);

%Stress recovery
[elementStressExtrapolated,elementStressAtGaussPoints]=stressRecext(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
[elementStressAtNodes]=stressRecovery_2(elementType,elements,nodes,constitutiveMatrix,displacementsVector);
% las tensiones extrapoladas son mas parecidas a las de NX

for iEle = 1:nElements
    elementNodes = nodes(elements(iEle,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = gaussPointsLocation(ipg,1);
        eta = gaussPointsLocation(ipg,2);
        
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
        % eleStress =  squeeze(stress(iEle,ipg,:));            % tensiones "directas"::: em los nodos
        eleStress =  squeeze(elementStressAtNodes(iEle,ipg,:));            % tensiones "directas"::: em los nodos

%         starStress = ( N * avgStress(elements(iEle,:),:) )'; % tensiones mejoradas::: en el centro tenes superconvergente del Q4.. 
%                                 promedias los centro de cada elemento a cada nodo y con eso tenes las tensiones mejoradas
        starStress = squeeze(elementStressExtrapolated(iEle,ipg,:)); % tensiones mejoradas        
        e2_el(iEle) = e2_el(iEle) + (starStress - eleStress)' * ... 
                    invC * (starStress - eleStress) * gaussPointsWeight(ipg) * det(jacobian);
                
        U2_el(iEle) = U2_el(iEle) + eleStress' * invC * eleStress * ...
                      gaussPointsWeight(ipg) * det(jacobian);
    end
    
    eta_el(iEle) = sqrt( e2_el(iEle) / (e2_el(iEle) + U2_el(iEle)) );
    
end

etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) );


%% Configuracion deformada
% D = (reshape(D,dofPerNode,[]))';
D = displacementsMatrix;
nodePosition = nodes + D(:,1:2);


%Graficacion
limites = [2.5E-1,3.5E-1];

labelNames = cellstr(num2str((1:size(nodes,1))'));

figure(1)
% bandPlot(elements,nodePosition,stressCentro(:,:,2),limites,'k');
bandPlot(elements,nodePosition,elementStressExtrapolated(:,:,2),limites,'k');
title('Tensiones en el centro de cada elemento.')
axis square

figure(2)
% bandPlot(elements,nodePosition,avgStress(:,2),limites,'k',[],'interp');
bandPlot(elements,nodePosition,avgStress(:,2),limites,'k',[],'interp');
title('Tensiones promediadas.')
axis square

figure(3)
bandPlot(elements,nodePosition,eta_el,[],'k',[],'flat');
title(' $\eta$ (Error  ZZ)','Interpreter','latex')
axis square

figure(4)
meshPlot(elements,nodePosition,'k');
text(nodePosition(:,1),nodePosition(:,2), labelNames)
title(' $\eta$ (Error  ZZ)','Interpreter','latex')
axis square