function [gaussPointsLocation,gaussPointsWeight] = getGaussPoints(elementShape,nGaussPoints)
% Gauss Point Location and Weight Calculator
% 
% [gaussPointsLocation,gaussPointsWeight] = gaussPoints(dimensions,nGaussPoints)
%
%  Gauss Points order
%
%  7---8---9        6
%  |       |        | \              
%  4   5   6        4  5   
%  |       |        |    \     
%  1---2---3        1--2--3

% gaussPointsLocation:  Gauss points locations in natural space cartesian coordinates
% gaussPointsWeight:    Gauss points weights
%
% elementShape:         Element Shape 'Line' 'Triangular' 'Quadrilateral' 'Tetrahedral' 'Hexahedral'
% nGaussPoints:         Number of Gauss Point of the scheme

%% Variables definition
gaussPointsWeight=zeros(nGaussPoints,1);
switch elementShape
    case {'Line'}
        gaussPointsLocation=zeros(nGaussPoints,1); n1DGaussPoints=nGaussPoints;
        numberOfGaussPointsAvailability=sum((n1DGaussPoints==[1 2 3 4 5]));
    case {'Triangular'}
        gaussPointsLocation=zeros(nGaussPoints,2);
        numberOfGaussPointsAvailability=sum((nGaussPoints==[1 3 4]));
    case {'Quadrilateral'}
        gaussPointsLocation=zeros(nGaussPoints,2); n1DGaussPoints=nGaussPoints^(1/2);
        numberOfGaussPointsAvailability=sum((n1DGaussPoints==[1 2 3 4 5]));
   case {'Tetrahedral'}
        gaussPointsLocation=zeros(nGaussPoints,3);
        numberOfGaussPointsAvailability=sum((nGaussPoints==[1 4 5]));
    case {'Hexahedral'}
        gaussPointsLocation=zeros(nGaussPoints,3); n1DGaussPoints=nGaussPoints^(1/3);
        numberOfGaussPointsAvailability=sum((n1DGaussPoints==[1 2 3 4 5]));
   otherwise
        disp('Element not recognized') %Controls if element type is available
        gaussPointsLocation=NaN;
        gaussPointsWeight=NaN;
        return
end

% Controls is number of points is available - Could not be available or not possible
if numberOfGaussPointsAvailability==0
    disp('Number of Gauss points not recognized')
    gaussPointsLocation=NaN;
    gaussPointsWeight=NaN;
    return
end


%Point location and weights calculations
switch elementShape
    case {'Line','Quadrilateral','Hexahedral'}
        switch n1DGaussPoints
            case 1
                linearGaussPointsWeight  = 2;
                linearGaussPointsLocation = 0;
            case 2
                linearGaussPointsWeight  = [1 1];
                a  = sqrt(3)/3;
                linearGaussPointsLocation = [-a a];
            case 3
                linearGaussPointsWeight  = [5/9 8/9 5/9];
                a  = sqrt(3/5);
                linearGaussPointsLocation = [-a 0 a];
            case 4
                a  = sqrt((3 - 2*sqrt(6/5))/7);
                b  = sqrt((3 + 2*sqrt(6/5))/7);
                linearGaussPointsLocation = [-b -a a b];
                wa = (18 + sqrt(30))/36;
                wb = (18 - sqrt(30))/36;
                linearGaussPointsWeight  = [wb wa wa wb];
            case 5
                a  = 1/3*sqrt(5 - 2*sqrt(10/7));
                b  = 1/3*sqrt(5 + 2*sqrt(10/7));
                linearGaussPointsLocation = [-b -a 0 a b];
                wa = (322 + 13*sqrt(70))/900;
                wb = (322 - 13*sqrt(70))/900;
                linearGaussPointsWeight  = [wb wa 128/225 wa wb];
        end
        switch elementShape
            case {'Line'}
                gaussPointsWeight=linearGaussPointsWeight;
                gaussPointsLocation=linearGaussPointsLocation;
            case {'Quadrilateral'}
                for i1DGaussPoint = 1:n1DGaussPoints
                        gaussPointsWeight((1:n1DGaussPoints)+(i1DGaussPoint-1)*n1DGaussPoints) = linearGaussPointsWeight(:)*linearGaussPointsWeight(i1DGaussPoint);
                        gaussPointsLocation((1:n1DGaussPoints)+(i1DGaussPoint-1)*n1DGaussPoints,1) = linearGaussPointsLocation(1:n1DGaussPoints);
                        gaussPointsLocation((1:n1DGaussPoints)+(i1DGaussPoint-1)*n1DGaussPoints,2) = linearGaussPointsLocation(i1DGaussPoint)*ones(n1DGaussPoints,1);
                end
            case {'Hexahedral'}
                for i1DGaussPoint = 1:n1DGaussPoints
                    for j1DGaussPoint = 1:n1DGaussPoints
                        gaussPointsWeight((1:n1DGaussPoints)+((j1DGaussPoint-1)+(i1DGaussPoint-1)*n1DGaussPoints)*n1DGaussPoints) = ...
                            linearGaussPointsWeight(:)*linearGaussPointsWeight(i1DGaussPoint)*linearGaussPointsWeight(j1DGaussPoint);
                        gaussPointsLocation((1:n1DGaussPoints)+((j1DGaussPoint-1)+(i1DGaussPoint-1)*n1DGaussPoints)*n1DGaussPoints,1) = ...
                            linearGaussPointsLocation(1:n1DGaussPoints);
                        gaussPointsLocation((1:n1DGaussPoints)+((j1DGaussPoint-1)+(i1DGaussPoint-1)*n1DGaussPoints)*n1DGaussPoints,2) = ...
                            linearGaussPointsLocation(j1DGaussPoint)*ones(n1DGaussPoints,1);
                        gaussPointsLocation((1:n1DGaussPoints)+((j1DGaussPoint-1)+(i1DGaussPoint-1)*n1DGaussPoints)*n1DGaussPoints,3) = ...
                            linearGaussPointsLocation(i1DGaussPoint)*ones(n1DGaussPoints,1);
                    end
                end
        end
    case {'Triangular'}
        switch nGaussPoints
            case 1
                gaussPointsWeight=1;
                gaussPointsLocation=[1 1]/3/2;
            case 3
                gaussPointsWeight=[1 1 1]'/3/2;
                gaussPointsLocation=[1 0 1;0 1 1]'/2;
            case 4
                gaussPointsWeight=[25 25 25 -27]'/48/2;
                gaussPointsLocation=[3/5 1/5 1/5 1/3; 1/3 1/5 3/5 1/3]';
        end
    case {'Tetrahedral'}
        switch nGaussPoints
            case 1
                gaussPointsWeight=1;
                gaussPointsLocation=[1 1 1]/4;
            case 4
                gaussPointsWeight=[1 1 1 1]'/4;
                a=(5+3*sqrt(5))/20; b=(5-sqrt(5))/20; 
                gaussPointsLocation=[a b b b; b b b a; b b a b]'/4;
            case 5
                gaussPointsWeight=[9 9 9 9 -16]'/20;
                gaussPointsLocation=[1/2 1/6 1/6 1/6 1/4; 1/6 1/6 1/6 1/2 1/4; 1/6 1/6 1/2 1/6 1/4]';
        end
end
