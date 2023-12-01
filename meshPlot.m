function meshPlot(elementNodesArray,nodesPositionArray,color,numbering)
% Mesh plotter
% 
% meshPlot(elementNodesArray,nodesPositionArray,color,numbering)
% 
% Nodal order in elementNodesArray:
%  Q4 Q8 Q9        CST LST
%
%  4---7---3        3
%  |       |        | \              
%  8   9   6        6  5   
%  |       |        |    \     
%  1---5---2        1--4--2
%  
% elementNodesArray:    element conectivity matrix
% nodesPositionArray:   nodal position in cartesian coordinates
% color:                element borders color
% numbering:            'Yes' 'No'
%

%% Input control
if (nargin < 4) || isempty(numbering)
    numbering='No';
end

%% Definition

nNodesPerElements = size(elementNodesArray,2);
nElements = size(elementNodesArray,1);

% Steering vector reorder definition
switch nNodesPerElements
    case {6}
        nodeReorderVector = [1 4 2 5 3 6];
    case {8,9}
        nodeReorderVector = [1 5 2 6 3 7 4 8];
    otherwise
        nodeReorderVector = 1:nNodesPerElements;
end

%% Object definition
h1 = patch('Faces',elementNodesArray(:,nodeReorderVector),'Vertices',nodesPositionArray);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
daspect([1 1 1]); hold on;

%% Entity numbering
if strcmp(numbering,'Yes')
    xCoordinateElementNodes=zeros(1,nNodesPerElements);
    yCoordinateElementNodes=zeros(1,nNodesPerElements);
    for iElements = 1:nElements;
        % Node numbering
        for iNodesPerElements = 1:nNodesPerElements;
            xCoordinateElementNodes(iNodesPerElements) = nodesPositionArray(elementNodesArray(iElements,iNodesPerElements),1);
            yCoordinateElementNodes(iNodesPerElements) = nodesPositionArray(elementNodesArray(iElements,iNodesPerElements),2);
            text(xCoordinateElementNodes(iNodesPerElements),yCoordinateElementNodes(iNodesPerElements),num2str(elementNodesArray(iElements,iNodesPerElements)),'VerticalAlignment','bottom','Color','r','FontSize',8);
        end
        plot(xCoordinateElementNodes,yCoordinateElementNodes,'o','MarkerFaceColor','r','MarkerEdgeColor','b','MarkerSize',4);
        % Element numbering
        xCoordinateElementCenter = mean(nodesPositionArray(elementNodesArray(iElements,:),1));
        yCoordinateElementCenter = mean(nodesPositionArray(elementNodesArray(iElements,:),2));
        text(xCoordinateElementCenter,yCoordinateElementCenter,num2str(iElements),'EdgeColor','k','Color','k','FontSize',8);
    end
end
