%% mallador

% cilindros zunchado
[elements1,nodes1,vertexNodes1,sideNodes]=quadrilateralDomainMeshGenerator_catedra(elementType,'Straight',b-a,h,nElementsR,nElementsZ,0,distorsion);
nodes1(:,1) = nodes1(:,1) + a;
[elements2,nodes2,vertexNodes2,sideNodes2] = quadrilateralDomainMeshGenerator(elementType,'Straight',c-b,h,nElementsR,nElementsZ,0,distorsion);
nNodes1 = size(nodes1,1); nElements1 = size(elements1,1); nNodes2 = size(nodes2,1);
vertexNodes2 = vertexNodes2 + nNodes1; 
sideNodes2 = sideNodes2 + size(nodes1,1);
nodes2(:,1) = nodes2(:,1) + b - interferencia ;
msh.elements = [elements1;elements2+size(nodes1,1)];
msh.nodes = [nodes1;nodes2];
msh.sideNodes = sideNodes;
msh.vertexNodes = vertexNodes1;

% esfera superior
meshAngle = 90;
meshInnerRadious = a;
meshLength = meshInnerRadious*meshAngle*pi/180;

[esf.elements,esf.nodes,esf.vertexNodes,esf.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,eEsf,50,nElementsR,meshAngle,0);

esf.nodes(:,2) = esf.nodes(:,2) + h;
meshPlot(esf.elements,esf.nodes,'r','No')

% merge mallas
tol = 0.001;
[msh] = mergeElAr(msh,esf,tol);

% Mesh plot
figure; meshPlot(msh.elements,msh.nodes,'b','No');