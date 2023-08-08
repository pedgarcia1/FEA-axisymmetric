function [stiffnessMatrix]=assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix)
% Finite element matrix assembler
% 
% assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray)
%
% stiffnessMatrix:      Assembled stiffness matrix
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
        nGaussPoints = 1;
        elementShape='Triangular';
    case 'LST'
        nDimensions = 2;                                %Problem dimension
        nGaussPoints = 3;
        elementShape='Triangular';
    case 'Q4'
        nDimensions = 2;                                %Problem dimension
        nGaussPoints = 4;
        elementShape='Quadrilateral';
    case 'Q8'
        nDimensions = 2;                                %Problem dimension
        nGaussPoints = 4;
        elementShape='Quadrilateral';
    case 'Q9'
        nDimensions = 2;                                %Problem dimension
        nGaussPoints = 9;
        elementShape='Quadrilateral';
end

nElements=size(elementNodesArray,1);            %Number of elements

nNodes=size(nodesPositionArray,1);              %Number of nodes
nElementalNodes = size(elementNodesArray,2);    %Number of node in each element
nElementalDof = nDimensions*nElementalNodes;    %Number of elemental Dofs

nTotalDof = nDimensions*nNodes;                 %Number of node in each element

%% Integration
[gaussPointsLocation,gaussPointsWeight] = getGaussPoints(elementShape,nGaussPoints);
gaussPointsLocation = gaussPointsLocation([1 2 4 3],:);

%% Stiffness matrix assembly
stiffnessMatrix = zeros(nTotalDof);

% Shape functions derivatives at every Gauss point in natural coordinates
shapeFunctionsDerivatives = getShapeFunctionsDerivatives(gaussPointsLocation,elementType);
[Ni,~] = shapefuns(gaussPointsLocation,elementType);

% q8 casero
% fN = @(ej,z) [- (ej^2*z)/4 + ej^2/4 - (ej*z^2)/4 + (ej*z)/4 + z^2/4 - 1/4, - (ej^2*z)/4 + ej^2/4 + (ej*z^2)/4 - (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/4 + ej^2/4 + (ej*z^2)/4 + (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/4 + ej^2/4 - (ej*z^2)/4 - (ej*z)/4 + z^2/4 - 1/4, (ej^2*z)/2 - z/2 - ej^2/2 + 1/2, ej/2 - (ej*z^2)/2 - z^2/2 + 1/2, z/2 - (ej^2*z)/2 - ej^2/2 + 1/2, (ej*z^2)/2 - ej/2 - z^2/2 + 1/2];
% dfN = @(ej,z) [ ej/2 + z/4 - (ej*z)/2 - z^2/4,  ej/2 - z/4 - (ej*z)/2 + z^2/4,  ej/2 + z/4 + (ej*z)/2 + z^2/4,  ej/2 - z/4 + (ej*z)/2 - z^2/4,    ej*z - ej, 1/2 - z^2/2,  - ej - ej*z, z^2/2 - 1/2;...
%               ej/4 + z/2 - (ej*z)/2 - ej^2/4, z/2 - ej/4 + (ej*z)/2 - ej^2/4, ej/4 + z/2 + (ej*z)/2 + ej^2/4, z/2 - ej/4 - (ej*z)/2 + ej^2/4, ej^2/2 - 1/2,  - z - ej*z, 1/2 - ej^2/2,    ej*z - z];
 
% for iGaussPoint = 1:nGaussPoints
%     Ni(:,:,iGaussPoint) = fN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
%     shapeFunctionsDerivatives(:,:,iGaussPoint) = dfN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
% end

% q4 casero 
% fN = @(ej,z) [(ej*z)/4 - z/4 - ej/4 + 1/4, ej/4 - z/4 - (ej*z)/4 + 1/4, ej/4 + z/4 + (ej*z)/4 + 1/4, z/4 - ej/4 - (ej*z)/4 + 1/4]; 
% dfN = @(ej,z) [ z/4 - 1/4,    1/4 - z/4,  z/4 + 1/4, - z/4 - 1/4 ; ej/4 - 1/4, - ej/4 - 1/4, ej/4 + 1/4,  1/4 - ej/4];
% gaussPointsLocation = gaussPointsLocation([1 2 4 3],:);
% for iGaussPoint = 1:nGaussPoints
%     Ni(:,:,iGaussPoint) = fN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
%     shapeFunctionsDerivatives(:,:,iGaussPoint) = dfN(gaussPointsLocation(iGaussPoint,1),gaussPointsLocation(iGaussPoint,2));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones de forma respecto de gaussPointLocation(iGaussPoint,1), gaussPointLocation(iGaussPoint,2)
%         N = shapefuns([gaussPointLocation(iGaussPoint,1) gaussPointLocation(iGaussPoint,2)],elementType);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iElement = 1:nElements
    elementalStiffnessMatrix = zeros(nElementalDof);
    
    elementalDofs = convertNode2Dof(elementNodesArray(iElement,:),nDimensions);
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:); 
    
    for iGaussPoint = 1:nGaussPoints
        
        r = Ni(:,:,iGaussPoint)*elementalNodesPosition(:,1);

        % Jacobian calculation at Gauss point
        jacobian = shapeFunctionsDerivatives(:,:,iGaussPoint)*elementalNodesPosition;
        jacobianDeterminant=det(jacobian);
        
        % Shape functions derivatives in structural coordinates
        structuralShapeFunctionsDerivatives = jacobian\shapeFunctionsDerivatives(:,:,iGaussPoint);
        
        B = zeros(size(constitutiveMatrix,1),nElementalDof);
        B(1,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(1,:);
        B(2,1:2:nElementalDof-1) =  Ni(:,:,iGaussPoint)/r; % N(xi,eta)_i/r(xi,eta) para cada pto de gauss
        B(3,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(2,:);
        B(4,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(2,:);
        B(4,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(1,:);
        
        elementalStiffnessMatrix = elementalStiffnessMatrix + 2*pi*B'*constitutiveMatrix*B*gaussPointsWeight(iGaussPoint)*jacobianDeterminant*r;
    end
    
    stiffnessMatrix(elementalDofs,elementalDofs) = stiffnessMatrix(elementalDofs,elementalDofs) + elementalStiffnessMatrix;
    
end

if any(diag(stiffnessMatrix) < 0)
    warning('elementos negativos en la diagonal. ')
end
end