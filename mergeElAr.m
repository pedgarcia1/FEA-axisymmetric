function [msh]=mergeElAr(stra,strb,tol)
nodar1=stra.nodesPositionArray;
elar=stra.elementNodesArray;
nodar2=strb.nodesPositionArray;
elarb=strb.elementNodesArray;
% sidea=stra.sideNodes;
% sideb=strb.sideNodes;


elnew=[elar;elarb+length(nodar1)];
nodnew=[nodar1;nodar2];
% isnoda=ismember(nodar1,nodar2,'rows');
% isnodb=ismember(nodar2,nodar1,'rows');
isnoda=ismembertol(nodar1,nodar2,tol,'ByRows',true);
isnodb=ismembertol(nodar2,nodar1,tol,'ByRows',true);

for i=1:length(sum(isnoda))
    nodac=find(isnoda);
    nodac=nodac(i);
    nodacb=find(isnodb);
    nodacb=nodacb(i);
    asd=ismember(elar,nodac);
    [fil,col ]=find(asd);
    elnew(fil,col)=nodacb+length(nodar1);
    
    aux = ismember(stra.sideNodes(4,:),nodac);
    [fil1,col1]=find(aux);
%     sidenewI(fil,col)=nodacb+length(nodar1);
    stra.sideNodes(4,col1) = 0;

    aux = ismember(stra.sideNodes(2,:),nodac);
    [fil2,col2]=find(aux);
%     sidenewD(fil,col)=nodacb+length(nodar1);
    stra.sideNodes(2,col2) = 0;
end

aux = strb.sideNodes(4,:);
aux(aux~=0) = aux(aux~=0)+length(nodar1);
log = stra.sideNodes(4,:)~=0;
sidenewI = [stra.sideNodes(4,log),aux];

aux = strb.sideNodes(2,:);
aux(aux~=0) = aux(aux~=0)+length(nodar1);
log = stra.sideNodes(2,:)~=0;
sidenewD = [stra.sideNodes(2,log),aux];

sidenew = zeros(4,max([length(sidenewI) length(sidenewD)]));
sidenew(4,1:length(sidenewI)) = sidenewI;
sidenew(2,1:length(sidenewD)) = sidenewD;

msh.elementNodesArray = elnew;
msh.nodesPositionArray = nodnew;
msh.sideNodes = sidenew;

% side nodes si noe s cero sumar length(nodar1) despues hacer otro for y si
% se repite el nodo (sin sumar) reemplazarlo por cero