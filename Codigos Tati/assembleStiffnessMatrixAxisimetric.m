function K=assembleStiffnessMatrixAxisimetric(Nodos,Elementos,C,eleType,eje)


% Gauss           
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
% Funciones de forma y derivadas
Ni = shapefuns   (upg,eleType);
dN = shapefunsder(upg,eleType);

nDofTot=max(size(Nodos)).*2;
nel=size(Elementos,1);
nNodEle=4;
nDofNod=2;
K = zeros(nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = Nodos(Elementos(iele,:),:);
    for ipg = 1:npg
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,ipg)*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,ipg);          % dNxy = inv(jac)*dN(:,:,ipg)
        

        

        if strcmp(eje,'y')
            r = Ni(:,:,ipg)*nodesEle(:,2);
        elseif strcmp(eje,'x')
            r = Ni(:,:,ipg)*nodesEle(:,1);
        end

        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = Ni(:,:,ipg)/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);        
        
        Ke = Ke + B'*C*B*r*det(jac)*wpg(ipg);
    end
    eleDofs = node2dof(Elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end