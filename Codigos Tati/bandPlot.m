function bandPlot(elementNodesArray,nodesPositionArray,variable,bandplotLimits,nBands,lineColor)
% Bandplot plotter
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
% variable:             Variable field to plot at nodes (Elements, Points, Values)
% bandplotLimits         Limits for the band plot
% lineColor:            Color for element borders
% nBands:               Number of band to divide bandplot        

%%  Definitions
error(nargchk(3, 6, nargin))
if (nargin < 4) || isempty(bandplotLimits)
    bandplotLimits = [min(min(variable)) max(max(variable))];
end

if (nargin < 5) || isempty(lineColor)
    lineColor = 'k';
end

tol = 1E-5;
nMinBandPlotDivisions = 3;
if (nargin < 6) || isempty(nBands)
    if diff(bandplotLimits) < bandplotLimits(2)*tol
        nBands = nMinBandPlotDivisions;
    else
        nBands = 10;
    end
end

if diff(bandplotLimits) < bandplotLimits(2)*tol
    bandplotLimits = bandplotLimits + [-1 1]*bandplotLimits(2)*tol;
end

if nBands < nMinBandPlotDivisions
    nBands = nMinBandPlotDivisions;
end

nElementNodes = size(elementNodesArray,2);
nElement = size(elementNodesArray,1);

switch nElementNodes
    case {6}
        reorderNodeNumbers = [1 4 2 5 3 6];
    case {8,9}
        reorderNodeNumbers = [1 5 2 6 3 7 4 8];
    otherwise
        reorderNodeNumbers = 1:nElementNodes;
end

% Patch creation
for iElement = 1:nElement
    elementNodes = elementNodesArray(iElement,:);
    h = patch('Faces',reorderNodeNumbers,'Vertices',nodesPositionArray(elementNodes,:),'FaceVertexCData',variable(iElement,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'CDataMapping','scaled');
end

% Scaling
colormap(jet(nBands))
caxis(bandplotLimits)

if nBands <= 20
    nTicks = nBands;
else
    nTicks = 20;
end

ticks = bandplotLimits(1):((diff(bandplotLimits))/nTicks):bandplotLimits(2);
tickLabels = cell(size(ticks));

for iTick = 1:length(ticks)
    tickLabels{iTick} = sprintf('%6.5E',ticks(iTick));
end

colorbar('YTick',ticks,'YTickLabel',tickLabels);
set(gca,'XTick',[],'YTick',[])
daspect([1 1 1])



