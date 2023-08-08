function [shapeFunctionsDerivatives] = getShapeFunctionsDerivatives(evaluationPointsArray,elementType)
% Evaluates derived shape functions at given points 
%
% [shapeFunctionsDerivatives] = getShapeFunctionsDerivatives(evaluationPointsArray,elementType)
%
% shapeFunctionsDerivatives: (Derivative,Function,Point)
%
% evaluationPointsArray:    Array with points where to evaluate shape functions derivatives (Coordinates,Points) 
% elementType:              Type of element 'CST' LST' 'Q4' 'Q8' 'Q9' 'AHMAD4' 'AHMAD8' 'AHMAD9'
%

%% Definitions
nPoints=size(evaluationPointsArray,1);
nDimensions=size(evaluationPointsArray,2);

% Natural coordinates
ksi=evaluationPointsArray(:,1);
eta=evaluationPointsArray(:,2);

switch elementType
    case {'Q4', 'AHMAD4'}
        shapeFunctionsDerivatives = zeros(2,4,nPoints);
        for iPoint = 1:nPoints
             shapeFunctionsDerivatives(:,:,iPoint) = [                  % ksi derivatives
                -0.25*(1 - eta(iPoint)),  0.25*(1 - eta(iPoint)), 0.25*(1 + eta(iPoint)), -0.25*(1 + eta(iPoint))
                                                                        % eta derivatives
                -0.25*(1 - ksi(iPoint)), -0.25*(1 + ksi(iPoint)), 0.25*(1 + ksi(iPoint)),  0.25*(1 - ksi(iPoint)) ];
        end
        
    case {'Q9', 'AHMAD9'}
          shapeFunctionsDerivatives = zeros(2,9,nPoints);
          for iPoint = 1:nPoints
            shapeFunctionsDerivatives(:,:,iPoint) = [                   % ksi derivatives
                0.25*eta(iPoint)*(-1+eta(iPoint))*(2*ksi(iPoint)-1),      0.25*eta(iPoint)*(-1+eta(iPoint))*(2*ksi(iPoint)+1),       0.25*eta(iPoint)*(1+eta(iPoint))*(2*ksi(iPoint)+1),...
                0.25*eta(iPoint)*( 1+eta(iPoint))*(2*ksi(iPoint)-1),                -ksi(iPoint)*eta(iPoint)*(-1+eta(iPoint)),  -1/2*(-1+eta(iPoint))*(1+eta(iPoint))*(2*ksi(iPoint)+1),...
                           -ksi(iPoint)*eta(iPoint)*(1+eta(iPoint)),  -1/2*(-1+eta(iPoint))*(1+eta(iPoint))*(2*ksi(iPoint)-1),           2*ksi(iPoint)*(-1+eta(iPoint))*(1+eta(iPoint))
                                                                        % eta derivatives
                   0.25*ksi(iPoint)*(-1+2*eta(iPoint))*(ksi(iPoint)-1),      0.25*ksi(iPoint)*(-1+2*eta(iPoint))*(1+ksi(iPoint)),       0.25*ksi(iPoint)*(2*eta(iPoint)+1)*(1+ksi(iPoint)),...
                    0.25*ksi(iPoint)*(2*eta(iPoint)+1)*(ksi(iPoint)-1),  -0.5*(ksi(iPoint)-1)*(1+ksi(iPoint))*(-1+2*eta(iPoint)),                 -ksi(iPoint)*eta(iPoint)*(1+ksi(iPoint)),...
                -0.5*(ksi(iPoint)-1)*(1+ksi(iPoint))*(2*eta(iPoint)+1),                 -ksi(iPoint)*eta(iPoint)*(ksi(iPoint)-1),           2*(ksi(iPoint)-1)*(1+ksi(iPoint))*eta(iPoint) ];
        end

    case {'Q8', 'AHMAD8'}
        shapeFunctionsDerivatives = zeros(2,8,nPoints);
        for iPoint = 1:nPoints
            shapeFunctionsDerivatives(:,:,iPoint) = [                   % ksi derivatives
                -0.25*(-1+eta(iPoint))*(eta(iPoint)+2*ksi(iPoint)),  -0.25*(-1+eta(iPoint))*(-eta(iPoint)+2*ksi(iPoint)),    0.25*(1+eta(iPoint))*(eta(iPoint)+2*ksi(iPoint)),   0.25*(1+eta(iPoint))*(-eta(iPoint)+2*ksi(iPoint)),...
                              ksi(iPoint)*(-1+eta(iPoint)),        -0.5*(-1+eta(iPoint))*(1+eta(iPoint)),                -ksi(iPoint)*(1+eta(iPoint)),        0.5*(-1+eta(iPoint))*(1+eta(iPoint))
                                                                        % eta derivatives
                -0.25*(-1+ksi(iPoint))*(ksi(iPoint)+2*eta(iPoint)),   -0.25*(1+ksi(iPoint))*(ksi(iPoint)-2*eta(iPoint)),    0.25*(1+ksi(iPoint))*(ksi(iPoint)+2*eta(iPoint)),   0.25*(-1+ksi(iPoint))*(ksi(iPoint)-2*eta(iPoint)),...
                      0.5*(-1+ksi(iPoint))*(1+ksi(iPoint)),                -(1+ksi(iPoint))*eta(iPoint),       -0.5*(-1+ksi(iPoint))*(1+ksi(iPoint)),               (-1+ksi(iPoint))*eta(iPoint) ];
        end
        
    case {'CST'}
        shapeFunctionsDerivatives = zeros(2,3,nPoints);
        for iPoint = 1:nPoints
            shapeFunctionsDerivatives(:,:,iPoint) = [-1, 1, 0           % ksi derivatives
                                                     -1, 0, 1];         % eta derivatives
        end
        
    case {'LST'}
        shapeFunctionsDerivatives = zeros(2,6,nPoints);
        for iPoint = 1:nPoints
            shapeFunctionsDerivatives(:,:,iPoint) = [
                 4*eta(iPoint)+4*ksi(iPoint)-3, 4*ksi(iPoint)- 1,               0, 4-8*ksi(iPoint)-4*eta(iPoint), 4*eta(iPoint),                -4*eta(iPoint)           % ksi derivatives
                 4*eta(iPoint)+4*ksi(iPoint)-3,                0, 4*eta(iPoint)-1,                -4*ksi(iPoint), 4*ksi(iPoint), 4-4*ksi(iPoint)-8*eta(iPoint)];         % eta derivatives
        end
        
end

