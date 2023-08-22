%% mallador
% elementType:          Type of element 'CST' LST' 'Q4' 'Q8' 'Q9'
% quadrilateralShape:   'Curved' 'Straight'
% meshLength:           Mesh length
% meshHeight:           Mesh height
% nNodesInLength:       Number of divisions in length 
% nNodesInHeight:       Number of divisions in height 
% meshAngle:            MeshAngle [degrees]
% distortion:           Regular distortion applied to the nodes position

%sideNodes              4lados x num nodos. lados son ABAJO DER ARRIBE IZQ
clear; close all;
elementType = 'Q4';
espesor = 50;
nElementsInLength = 8;
set(0,'DefaultFigureWindowStyle','docked')
% deltaEsp = 15;
nEleConc1 = 1;
nEleConc2 = 5;

lengthCilindro = 1100;

%% cilindro
[cil.elementNodesArray,cil.nodesPositionArray,cil.vertexNodes,cil.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Straight',espesor,lengthCilindro,nElementsInLength,100,0,0);

cil.nodesPositionArray(:,1) = cil.nodesPositionArray(:,1) + 300;


figure
meshPlot(cil.elementNodesArray,cil.nodesPositionArray,'b','No');
hold on

%% cilindro concentrador 
[con.elementNodesArray,con.nodesPositionArray,con.vertexNodes,con.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Straight',espesor + espesor/nElementsInLength*nEleConc1,1200-lengthCilindro,nElementsInLength + nEleConc1,15,0,0);

con.nodesPositionArray(:,1) = con.nodesPositionArray(:,1) + 300;
con.nodesPositionArray(:,2) = con.nodesPositionArray(:,2) + lengthCilindro;


% figure
meshPlot(con.elementNodesArray,con.nodesPositionArray,'b','No');

%% esfera sup
meshAngle = 90;
meshInnerRadious = 300;
meshLength = meshInnerRadious*meshAngle*pi/180;

[esf.elementNodesArray,esf.nodesPositionArray,esf.vertexNodes,esf.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,espesor,50,nElementsInLength,meshAngle,0);

esf.nodesPositionArray(:,2) = esf.nodesPositionArray(:,2) + 1200;

% figure
meshPlot(esf.elementNodesArray,esf.nodesPositionArray,'b','No');

%% conc 2 abajo
meshAngle = 45;
meshInnerRadious = 300;
meshLength = meshInnerRadious*meshAngle*pi/180;

[con2.elementNodesArray,con2.nodesPositionArray,con2.vertexNodes,con2.sideNodes]=quadrilateralDomainMeshGenerator_2(elementType,'Curved',meshLength,espesor + espesor/nElementsInLength*nEleConc2,20,nElementsInLength + nEleConc2,meshAngle,0);
% [esf.elementNodesArray,esf.nodesPositionArray,esf.vertexNodes,esf.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,espesor,20,nElementsInLength,meshAngle,0);

con2.nodesPositionArray(:,2) = con2.nodesPositionArray(:,2);

figure
meshPlot(con2.elementNodesArray,con2.nodesPositionArray,'b','No');

%% recta inf

[inf.elementNodesArray,inf.nodesPositionArray,inf.vertexNodes,inf.sideNodes] = quadrilateralDomainMeshGenerator(elementType,'Straight',sqrt(2)*300,espesor,30,nElementsInLength,0);

theta = -pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

inf.nodesPositionArray_rotated = inf.nodesPositionArray * R;
inf.nodesPositionArray = inf.nodesPositionArray_rotated;

outLado = inf.sideNodes(2,:);
inLado = inf.sideNodes(4,:);

figure
meshPlot(inf.elementNodesArray,inf.nodesPositionArray,'b','No');

% p1 = nodesPositionArray(nodosLado(1),:);
% aux = nodosLado(nodosLado~=0);
% p2 = nodesPositionArray(aux(end),:);
% angle = atan2(p2(2)-p2(1), p2(1)-p1(1)) - pi/2; % da 45 boludazo

%% radio
meshAngle = 45;
meshInnerRadious = 300;
meshLength = meshInnerRadious*meshAngle*pi/180;

[rad.elementNodesArray,rad.nodesPositionArray,rad.vertexNodes,rad.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,espesor + espesor/nElementsInLength*nEleConc2,20,nElementsInLength + nEleConc2,meshAngle,0);
% [esf.elementNodesArray,esf.nodesPositionArray,esf.vertexNodes,esf.sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Curved',meshLength,espesor,20,nElementsInLength,meshAngle,0);

rad.nodesPositionArray(:,2) = -rad.nodesPositionArray(:,2);

figure
meshPlot(rad.elementNodesArray,rad.nodesPositionArray,'b','No');


%%
tol = 0.001;
[msh] = mergeElAr(esf,con,tol);
[msh] = mergeElAr(msh,cil,tol);
[msh] = mergeElAr(msh,con2,tol);

figure
meshPlot(msh.elementNodesArray,msh.nodesPositionArray,'b','No')
