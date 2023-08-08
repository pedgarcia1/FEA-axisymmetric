function [elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator_2(elementType,quadrilateralShape,meshLength,meshHeight,nElementsInLength,nElementsInHeight,meshAngle,distortion)
% Mesh Generator for a Quadrilateral Domain
% 
% [nodesPositionArray,elementNodesArray,vertexNode,sideNodes]=quadrilateralDomainMeshGenerator(elementType,quadrilateralShape,meshLength,meshHeight,nElementsInLength,nElementsInHeight,meshAngle)
%
% nodesPositionArray:   Nodal position in cartesian coordinates
% elementNodesArray:    Element conectivity matrix
% vertexNode:           Node number of the mesh vertex points
% sideNodes:            Nodes numbers of the mesh sides
%
% elementType:          Type of element 'CST' LST' 'Q4' 'Q8' 'Q9'
% quadrilateralShape:   'Curved' 'Straight'
% meshLength:           mMesh length
% meshHeight:           Mesh height
% nNodesInLength:       Number of divisions in length 
% nNodesInHeight:       Number of divisions in height 
% meshAngle:            MeshAngle [degrees]
% distortion:           Regular distortion applied to the nodes position
%

%% Input parameter check

if nargin<=5
    nElementsInLength=1; nElementsInHeight=1;  meshAngle=0; distortion=0;
elseif nargin<7
    meshAngle=0; distortion=0;
elseif nargin<8
    distortion=0;
end

if meshAngle==0
    quadrilateralShape='Straight';
end

%% Variables definition
switch elementType
    case 'CST'
        nNodEle = 3;                                        % Nodes per element
        nNodesInLength=nElementsInLength+1;   nNodesInHeight=nElementsInHeight+1;
    case 'LST'
        nNodEle = 6;                                        % Nodes per element
        nNodesInLength=2*nElementsInLength+1; nNodesInHeight=2*nElementsInHeight+1; 
    case 'Q4'
        nNodEle = 4;                                        % Nodes per element
        nNodesInLength=nElementsInLength+1;   nNodesInHeight=nElementsInHeight+1;
    case {'Q8','Q9'}
        nNodEle = 9;                                        % Nodes per element 
        nNodesInLength=2*nElementsInLength+1; nNodesInHeight=2*nElementsInHeight+1;
end


%% Node placement
switch quadrilateralShape
    case 'Straight'; % Nodal Discretization Rectangular Mesh
        [horizontalPosition,verticalPosition]=meshgrid((0:meshHeight/(nNodesInHeight-1):meshHeight),0:meshLength/(nNodesInLength-1):meshLength);
        %Mesh distortion
        horizontalPosition(2:nNodesInLength-1,2:nNodesInHeight-1)=horizontalPosition(2:nNodesInLength-1,2:nNodesInHeight-1)+distortion.*(rand(nNodesInLength-2,nNodesInHeight-2)-0.5);
        verticalPosition(2:nNodesInLength-1,2:nNodesInHeight-1)=verticalPosition(2:nNodesInLength-1,2:nNodesInHeight-1)+distortion.*(rand(nNodesInLength-2,nNodesInHeight-2)-0.5);
        %Coordinates asignment to nodes
        nodesPositionArray(:,1)=verticalPosition(:);
        nodesPositionArray(:,2)=horizontalPosition(:);
    case 'Curved'; % Nodal Discretization Circular Mesh
        meshInnerRadious=meshLength/meshAngle*180/pi;
        [radialPosition,angularPosition]=meshgrid(meshInnerRadious+(0:meshHeight/(nNodesInHeight-1):meshHeight),90:meshAngle/(nNodesInLength-1):(90+meshAngle));
        %Mesh distortion
        radialPosition(2:nNodesInLength-1,2:nNodesInHeight-1)=radialPosition(2:nNodesInLength-1,2:nNodesInHeight-1)+distortion.*(rand(nNodesInLength-2,nNodesInHeight-2)-0.5);
        angularPosition(2:nNodesInLength-1,2:nNodesInHeight-1)=angularPosition(2:nNodesInLength-1,2:nNodesInHeight-1)+distortion.*(rand(nNodesInLength-2,nNodesInHeight-2)-0.5);
        %Coordinates asignment to nodes
        nodesPositionArray(:,1)=radialPosition(:).*sind(angularPosition(:));
        nodesPositionArray(:,2)=radialPosition(:).*cosd(angularPosition(:));
end

%% Connectivity Matrix
% Number of elements calculation
switch elementType
    case {'CST','LST'}
        nElements=2*nElementsInLength*nElementsInHeight;
    case {'Q4','Q8','Q9'}
        nElements=nElementsInLength*nElementsInHeight;
end
elementNodesArray = zeros(nElements,nNodEle);

% First element structure calculation
switch elementType
    case 'CST'
        elementNodesArray(1,:)=[1     2                    nNodesInLength+1 ];
        elementNodesArray(2,:)=[2     nNodesInLength+2     nNodesInLength+1 ];
    case 'LST'
        elementNodesArray(1,:)=[1                   3                     2*nNodesInLength+1 ...
                                2                   nNodesInLength+2      nNodesInLength+1];
        elementNodesArray(2,:)=[3                   2*nNodesInLength+3    2*nNodesInLength+1 ...
                                nNodesInLength+3    2*nNodesInLength+2      nNodesInLength+2];
    case 'Q4'
        elementNodesArray(1,:)=[1     2                    nNodesInLength+2     nNodesInLength+1 ];
    case 'Q8'
        elementNodesArray(1,:)=[1     3                    2*nNodesInLength+3-nElementsInLength   2*nNodesInLength+1-nElementsInLength ...
                                2     nNodesInLength+2     2*nNodesInLength+2-nElementsInLength   nNodesInLength+1   ...
                                nNodesInLength+2]; % Central node stored to know which nodes to eliminate in nodesPositionArray
    case 'Q9'
        elementNodesArray(1,:)=[1     3                    2*nNodesInLength+3   2*nNodesInLength+1 ...
                                2     nNodesInLength+3     2*nNodesInLength+2   nNodesInLength+1   ...
                               nNodesInLength+2];
end
% Element displacement acording to element number
iElement=0;
for iElementsInHeight=1:nElementsInHeight
    for iElementsInLength=1:nElementsInLength
        switch elementType
            case 'CST'
                iElement=iElement+2;
                elementNodesArray(iElement-1,:)=elementNodesArray(1,:)+(iElementsInLength-1)+(iElementsInHeight-1)*nNodesInLength;
                elementNodesArray(iElement,:)=elementNodesArray(2,:)+(iElementsInLength-1)+(iElementsInHeight-1)*nNodesInLength;
            case 'LST'
                iElement=iElement+2;
                elementNodesArray(iElement-1,:)=elementNodesArray(1,:)+2*(iElementsInLength-1)+(iElementsInHeight-1)*2*nNodesInLength;
                elementNodesArray(iElement,:)=elementNodesArray(2,:)+2*(iElementsInLength-1)+(iElementsInHeight-1)*2*nNodesInLength;
            case 'Q4'
                iElement=iElement+1;
                elementNodesArray(iElement,:)=elementNodesArray(1,:)+(iElementsInLength-1)+(iElementsInHeight-1)*nNodesInLength;
            case 'Q8'
                iElement=iElement+1;
                elementNodesArray(iElement,[1 2 3 4 5 7])=elementNodesArray(1,[1 2 3 4 5 7])+2*(iElementsInLength-1)+(iElementsInHeight-1)*(2*nNodesInLength-nElementsInLength);
                elementNodesArray(iElement,[6 8])=elementNodesArray(1,[6 8])+(iElementsInLength-1)+(iElementsInHeight-1)*(2*nNodesInLength-nElementsInLength);
                elementNodesArray(iElement,9)=elementNodesArray(1,9)+2*(iElementsInLength-1)+(iElementsInHeight-1)*2*nNodesInLength;
            case 'Q9'
                iElement=iElement+1;
                elementNodesArray(iElement,:)=elementNodesArray(1,:)+2*(iElementsInLength-1)+(iElementsInHeight-1)*2*nNodesInLength;
        end
    end
end

if strcmp(elementType,'Q8');
    nodesPositionArray(elementNodesArray(:,9),:)=[];    %Unused nodes elimination
    elementNodesArray(:,9)=[] ;                         %Unused central node elimination
end

%% Boundaries listing

switch elementType
    case {'CST','LST','Q4','Q9'}
        %Vertex nodes
        vertexNodes(1)=1;                                    %South West
        vertexNodes(2)=nNodesInLength;                       %South East
        vertexNodes(3)=(nNodesInHeight-1)*nNodesInLength+1;  %North West
        vertexNodes(4)=nNodesInLength*nNodesInHeight;        %North East
        
        %Side nodes
        sideNodes=zeros(4,max(nNodesInLength,nNodesInHeight));
        sideNodes(1,1:nNodesInLength)=1:nNodesInLength;                                                        %South
        sideNodes(2,1:nNodesInHeight)=nNodesInLength:nNodesInLength:nNodesInLength*nNodesInHeight;             %East
        sideNodes(3,1:nNodesInLength)=((nNodesInHeight-1)*nNodesInLength+1):(nNodesInLength*nNodesInHeight);   %North
        sideNodes(4,1:nNodesInHeight)=1:nNodesInLength:((nNodesInHeight-1)*nNodesInLength+1);                  %West
    case {'Q8'}
        %Vertex nodes
        vertexNodes(1)=1;                                                            %South West
        vertexNodes(2)=nNodesInLength;                                               %South East
        vertexNodes(3)=nNodesInLength*(nNodesInHeight-1)+1-nElements;                %North West
        vertexNodes(4)=nNodesInLength*nNodesInHeight-nElements;                      %North East
        
        %Side nodes
        sideNodes=zeros(4,max(nNodesInLength,nNodesInHeight));
        sideNodes(1,1:nNodesInLength)=1:nNodesInLength;                                                                                                                     %South
        sideNodes(2,1:2:nNodesInHeight)=nNodesInLength:2*nNodesInLength-nElementsInLength:nNodesInLength*nNodesInHeight-nElements;                                          %East
        sideNodes(2,2:2:nNodesInHeight-1)=2*nNodesInLength-nElementsInLength:2*nNodesInLength-nElementsInLength:(2*nNodesInLength-nElementsInLength)*nElementsInHeight;     %East
        sideNodes(3,1:nNodesInLength)=(nNodesInLength*(nNodesInHeight-1)+1-nElements):nNodesInLength*nNodesInHeight-nElements;                                              %North
        sideNodes(4,1:2:nNodesInHeight)=1:2*nNodesInLength-nElementsInLength:nNodesInLength*(nNodesInHeight-1)+1-nElements;                                                 %West
        sideNodes(4,2:2:nNodesInHeight-1)=nNodesInLength+1:2*nNodesInLength-nElementsInLength:nNodesInLength+1+(2*nNodesInLength-nElementsInLength)*(nElementsInHeight-1);  %East
end

%% Control Plot

% figure
% meshPlot(elementNodesArray,nodesPositionArray,'b','No');

return