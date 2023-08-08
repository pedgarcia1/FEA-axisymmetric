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

%% Stiffness matrix assembly
stiffnessMatrix = zeros(nTotalDof);

% Shape functions derivatives at every Gauss point in natural coordinates
shapeFunctionsDerivatives = getShapeFunctionsDerivatives(gaussPointsLocation,elementType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones de forma respecto de gaussPointLocation(iGaussPoint,1), gaussPointLocation(iGaussPoint,2)
%         N = shapefuns([gaussPointLocation(iGaussPoint,1) gaussPointLocation(iGaussPoint,2)],elementType);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iElement = 1:nElements
    elementalStiffnessMatrix = zeros(nElementalDof);
    
    elementalDofs = convertNode2Dof(elementNodesArray(iElement,:),nDimensions);
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:);
    for iGaussPoint = 1:nGaussPoints
        
        % Jacobian calculation at Gauss point
        jacobian = shapeFunctionsDerivatives(:,:,iGaussPoint)*elementalNodesPosition;
        jacobianDeterminant=det(jacobian);
        
        % Shape functions derivatives in structural coordinates
        structuralShapeFunctionsDerivatives = jacobian\shapeFunctionsDerivatives(:,:,iGaussPoint);
        
        B = zeros(size(constitutiveMatrix,1),nElementalDof);
        B(1,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(1,:);
        B(2,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(2,:);
        B(3,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(2,:);
        B(3,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(1,:);
        
        elementalStiffnessMatrix = elementalStiffnessMatrix + B'*constitutiveMatrix*B*gaussPointsWeight(iGaussPoint)*jacobianDeterminant;
    end
    
    stiffnessMatrix(elementalDofs,elementalDofs) = stiffnessMatrix(elementalDofs,elementalDofs) + elementalStiffnessMatrix;
end

