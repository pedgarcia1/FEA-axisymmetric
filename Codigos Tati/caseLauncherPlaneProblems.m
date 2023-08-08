% Case launcher
clear; close all; 

%% Preprocess

elementType='Q9';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
loadCase='Constant Bending';         %'Uniform' 'Constant Bending' 'Variable Bending'

% Mesh generation
[elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator(elementType,'Straight',2,1,1,1,45,0);
    
nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

% Material properties
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,1,1/3);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(vertexNodes(1),1:2) = true;
boundaryConditionsArray(sideNodes(end,find(sideNodes(end,:))),1) = true;

% Load definition
distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction

% Load case load imposition
switch loadCase
    case 'Uniform'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),1) = 1;
                pointLoadsArray (vertexNodes(2),1) = 1/2;
                pointLoadsArray (vertexNodes(4),1) = 1/2;
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,1:end),1) = 1/3;
                pointLoadsArray (sideNodes(2,2:2:end),1) = 2/3;
                pointLoadsArray (vertexNodes(2),1) = 1/6;
                pointLoadsArray (vertexNodes(4),1) = 1/6;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/abs(sum(pointLoadsArray(:,1)));
    case 'Constant Bending'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),1) = -1:2/(size(sideNodes(2,:),2)-1):1;
                pointLoadsArray (vertexNodes(2),1) = -1/2;
                pointLoadsArray (vertexNodes(4),1) =  1/2;
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,:),1) = (-1:2/(size(sideNodes(2,:),2)-1):1)/3;
                pointLoadsArray (sideNodes(2,2:2:end),1) = 2*pointLoadsArray (sideNodes(2,2:2:end),1);
                pointLoadsArray (vertexNodes(2),1) = -1/6;
                pointLoadsArray (vertexNodes(4),1) =  1/6;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/sum(pointLoadsArray(:,1).*nodesPositionArray(:,2));
    case 'Variable Bending'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),2) = -((-1:2/(size(sideNodes(2,:),2)-1):1).*(-1:2/(size(sideNodes(2,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(4,:),2) = ((-1:2/(size(sideNodes(4,:),2)-1):1).*(-1:2/(size(sideNodes(4,:),2)-1):1)+1);
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                % Right side shear
                pointLoadsArray (sideNodes(2,:),2) = -((-1:2/(size(sideNodes(2,:),2)-1):1).*(-1:2/(size(sideNodes(2,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(2,2:2:end),2) = 2*pointLoadsArray (sideNodes(2,2:2:end),1);
                pointLoadsArray (vertexNodes(2),1) = 0;
                pointLoadsArray (vertexNodes(4),1) =  0;
                % Left side shear
                pointLoadsArray (sideNodes(4,:),2) = ((-1:2/(size(sideNodes(4,:),2)-1):1).*(-1:2/(size(sideNodes(4,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(4,2:2:end),2) = 2*pointLoadsArray (sideNodes(4,2:2:end),1);
                pointLoadsArray (vertexNodes(1),1) = 0;
                pointLoadsArray (vertexNodes(3),1) =  0;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/abs(sum(pointLoadsArray(sideNodes(4,:),2)))/10;
end



%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix);

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = stiffnessMatrix(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;

%% Postprocess
%Stress recovery
[elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);

% Deformed plot
meshPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)','b','Yes');

% Stresses plot
bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));