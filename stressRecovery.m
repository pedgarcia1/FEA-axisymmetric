function [elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,nodalDisplacements)
% Stress recovery from displacements
% 
% [elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacements)
%
% elementStressAtNodes:      Element stresses at element nodes
%
% elementType:          Type of element 'CST' LST' 'Q4' 'Q8' 'Q9'
% nodesPositionArray:   Nodal position in cartesian coordinates
% elementNodesArray:    Element conectivity matrix
% constitutiveMatrix:   Constitutive Matrix
%

%% Definitions
switch elementType
    case 'CST'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates=[0 0
                                 1 0
                                 0 1];
    case 'LST'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates=[0   0
                                 1   0
                                 0   1
                                 0.5 0
                                 0.5 0.5
                                 0   1];
    case 'Q4'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1];
    case 'Q8'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1
                                    0 -1
                                    1  0
                                    0  1
                                   -1  0];
    case 'Q9'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1
                                    0 -1
                                    1  0
                                    0  1
                                   -1  0
                                    0  0]*sqrt(3)/3;
end

nElements=size(elementNodesArray,1);            %Number of elements

nNodes=size(nodesPositionArray,1);              %Number of nodes
nElementalNodes = size(elementNodesArray,2);    %Number of node in each element
nElementalDof = nDimensions*nElementalNodes;    %Number of elemental Dofs

nTotalDof = nDimensions*nNodes;                 %Number of node in each element

nPoints = size(naturalNodalCoordinates,1);      %Number of points to get stress

elementStressAtNodes = zeros(nElements,nPoints,4);

%% Stress recovery

% Shape functions derivatives at every point in natural coordinates
% shapeFunctionsDerivatives = getShapeFunctionsDerivatives(naturalNodalCoordinates,elementType);

% [gaussPointsLocation,] = getGaussPoints('Quadrilateral',4);
% gaussPointsLocation = gaussPointsLocation([1 2 4 3],:);

fN = @(ej,z) [- (ej^2*z)/4 + ej^2/4 - (ej*z^2)/4 + (ej*z)/4 + z^2/4 - 1/4, - (ej^2*z)/4 + ej^2/4 + (ej*z^2)/4 - (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/4 + ej^2/4 + (ej*z^2)/4 + (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/4 + ej^2/4 - (ej*z^2)/4 - (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/2 - z/2 - ej^2/2 + 1/2, ej/2 - (ej*z^2)/2 - z^2/2 + 1/2, z/2 - (ej^2*z)/2 - ej^2/2 + 1/2, (ej*z^2)/2 - ej/2 - z^2/2 + 1/2];
dfN = @(ej,z) [ ej/2 + z/4 - (ej*z)/2 - z^2/4,  ej/2 - z/4 - (ej*z)/2 + z^2/4,  ej/2 + z/4 + (ej*z)/2 + z^2/4,  ej/2 - z/4 + (ej*z)/2 - z^2/4,    ej*z - ej, 1/2 - z^2/2,  - ej - ej*z, z^2/2 - 1/2;...
              ej/4 + z/2 - (ej*z)/2 - ej^2/4, z/2 - ej/4 + (ej*z)/2 - ej^2/4, ej/4 + z/2 + (ej*z)/2 + ej^2/4, z/2 - ej/4 - (ej*z)/2 + ej^2/4, ej^2/2 - 1/2,  - z - ej*z, 1/2 - ej^2/2,    ej*z - z];
 
% for iGaussPoint = 1:nPoints
%     
%     Ni(:,:,iGaussPoint) = fN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
%     shapeFunctionsDerivatives(:,:,iGaussPoint) = dfN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
% 
% end


for iElement = 1:nElements
    
    elementalDofs = convertNode2Dof(elementNodesArray(iElement,:),nDimensions);
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:);
    for iPoint = 1:nPoints
        
        r = fN(naturalNodalCoordinates(iPoint,1),naturalNodalCoordinates(iPoint,2))*elementalNodesPosition(:,1);

        % Jacobian calculation at Gauss point
%         jacobian = shapeFunctionsDerivatives(:,:,iPoint)*elementalNodesPosition;
        jacobian = dfN(naturalNodalCoordinates(iPoint,1),naturalNodalCoordinates(iPoint,2))*elementalNodesPosition;
%         jacobianDeterminant=det(jacobian);
        
        % Shape functions derivatives in structural coordinates
%         structuralShapeFunctionsDerivatives = jacobian\shapeFunctionsDerivatives(:,:,iPoint);
        structuralShapeFunctionsDerivatives = jacobian\dfN(naturalNodalCoordinates(iPoint,1),naturalNodalCoordinates(iPoint,2));
        
        B = zeros(size(constitutiveMatrix,1),nElementalDof);
        B(1,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(1,:);
        B(2,1:2:nElementalDof-1) =  fN(naturalNodalCoordinates(iPoint,1),naturalNodalCoordinates(iPoint,2))/r;
        B(3,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(2,:);
        B(4,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(2,:);
        B(4,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(1,:);
        
        elementStressAtNodes(iElement,iPoint,:) = constitutiveMatrix*B*nodalDisplacements(elementalDofs);
    end
end

