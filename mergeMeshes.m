function [nodes, elements] = mergeMeshes(nodes1,elements1,nodes2,elements2)


nNod1=size(nodes1,1);
elem2BU=elements2;
elements2=elements2+nNod1;
%% Buscamos que nodos estan ocupando el mismo lugar en ambas mallas


NodosCompartidos=[];
tolerancia=1e-2;
for i = 1:size(nodes1, 1)
    nodo_malla1 = nodes1(i, :);    
    for j = 1:size(nodes2, 1)
        nodo_malla2 = nodes2(j, :);
        
        distancia = norm(nodo_malla1 - nodo_malla2);
        
        if distancia < tolerancia
            NodosCompartidos = [NodosCompartidos; i, j];
        end
    end
end



for i = 1:size(NodosCompartidos, 1)
    mask = elem2BU == NodosCompartidos(i, 2);
    elements2(mask) = NodosCompartidos(i, 1);
end



nodes2(NodosCompartidos(:, 2),:)=-ones(size(NodosCompartidos,1),2);




for k = 1:size(nodes2, 1)
    if all(nodes2(k, :) == -1)
        continue;
    end
    cantidadDesplazada = sum(nodes2(1:k, 1) == -1);
    mask = elem2BU == k;
    elements2(mask) = elements2(mask) - cantidadDesplazada;
end

nodes2(NodosCompartidos(:, 2),:)=[];

nodes=[nodes1;nodes2];
elements=[elements1;elements2];

end