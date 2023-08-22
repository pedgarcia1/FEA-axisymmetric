function [msh]=mergeElAr(stra,strb,tol)

nodosRepetidosB = find(ismembertol(strb.nodes,stra.nodes,tol,'ByRows',true));

strb.elements = strb.elements + size(stra.nodes,1);

for i = nodosRepetidosB'
    [nodoEnA] = find(ismembertol(stra.nodes,strb.nodes(i,:),tol,'ByRows',true));
    strb.elements(find(ismember(strb.elements,i+size(stra.nodes,1)))) = nodoEnA;
end

strb.nodes(nodosRepetidosB,:) = [];
msh.nodes = [stra.nodes;strb.nodes];
msh.elements = [stra.elements;strb.elements];
