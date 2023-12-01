function [side] =side2element(Side,elementNodesArray,sideNodes)

switch Side
    case 'Up'
        Position=3;
    case 'Down'
        Position=1;
    case 'Left'
        Position=4;
    case 'Right'
        Position=2;
end

[row,~]=find(ismember(elementNodesArray,sideNodes(Position,:)));
side=unique(row)';
end

