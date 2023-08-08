function [elementStressAtNodes]=stressRecext_2(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,nodalDisplacements)
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

nDimensions = 2;                                %Problem dimension
% naturalNodalCoordinates=[-1  -1
%     1  -1
%     1   1
%     -1   1]/sqrt(3);
naturalNodalCoordinates = [-1 -1
    1 -1
    1  1
    -1  1
    0 -1
    1  0
    0  1
    -1  0];

nElements=size(elementNodesArray,1);            %Number of elements

nNodes=size(nodesPositionArray,1);              %Number of nodes
nElementalNodes = size(elementNodesArray,2);    %Number of node in each element
nElementalDof = nDimensions*nElementalNodes;    %Number of elemental Dofs

nTotalDof = nDimensions*nNodes;                 %Number of node in each element

nPoints = size(naturalNodalCoordinates,1);      %Number of points to get stress

elementStressAtNodes = zeros(nElements,nPoints,4);

%% Stress recovery

% Shape functions derivatives at every point in natural coordinates
shapeFunctionsDerivatives = getShapeFunctionsDerivatives(naturalNodalCoordinates,elementType);

[gaussPointsLocation,~] = getGaussPoints('Quadrilateral',4);
gaussPointsLocation = gaussPointsLocation([1 3 2 4],:);
[Ni,~] = shapefuns(gaussPointsLocation,elementType);

for iElement = 1:nElements
    
    elementalDofs = convertNode2Dof(elementNodesArray(iElement,:),nDimensions);
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:);
    for iPoint = 1:nPoints
        
         r = Ni(:,:,iPoint)*elementalNodesPosition(:,1);

        % Jacobian calculation at Gauss point
        jacobian = shapeFunctionsDerivatives(:,:,iPoint)*elementalNodesPosition;
        jacobianDeterminant=det(jacobian);
        
        % Shape functions derivatives in structural coordinates
        structuralShapeFunctionsDerivatives = jacobian\shapeFunctionsDerivatives(:,:,iPoint);
        
        B = zeros(size(constitutiveMatrix,1),nElementalDof);
        B(1,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(1,:);
        B(2,1:2:nElementalDof-1) =  Ni(:,:,iPoint)/r; % N(xi,eta)_i/r(xi,eta) para cada pto de gauss
        B(3,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(2,:);
        B(4,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(2,:);
        B(4,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(1,:);
        
        elementStressAtNodes(iElement,iPoint,:) = constitutiveMatrix*B*nodalDisplacements(elementalDofs);
    end
    
end

% Extrapolo a los Nodos

a = sqrt(3);
rsExt = a*[-1  -1
    1  -1
    1   1
    -1   1
    0  -1
    1   0
    0   1
    -1   0];
nPoints = size(rsExt,1); 
stressExtra = zeros(nElements,nPoints,4); 
% SON 4 PORQUE SON sigmar sigmatit sigmaz y tau_zr
for iele = 1:nElements
    for inode = 1:nPoints
        r = rsExt(inode,1);
        s = rsExt(inode,2);
        
        N = shapefuns([r s],'Q4');
        
        stressExtra(iele,inode,:) = N * squeeze(elementStressAtNodes(iele,:,:));
    end
end
end
