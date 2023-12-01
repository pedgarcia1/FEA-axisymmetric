function [pointLoadsArray] = DistributiveLoads(Sigma,elementType,elementNodesArray,nodesPositionArray,Side,sideNodes)

nDimensions=2;
nElements=size(elementNodesArray,1);            %Number of elements
nNodes=size(nodesPositionArray,1);              %Number of nodes
nElementalNodes = size(elementNodesArray,2);    %Number of node in each element
nElementalDof = nDimensions*nElementalNodes;    %Number of elemental Dofs
nTotalDof = nDimensions*nNodes;                 %Number of node in each element


%% Loads
pointLoadsArray = zeros(nNodes,nDimensions); 
pointLoadsArray2 = zeros(nNodes,nDimensions); 

elementVector=side2element(Side,elementNodesArray,sideNodes);
[PointIntegration,gaussPointsWeight] =getPointInteger(Side,elementType);

nPointInterger=size(PointIntegration,1);

for iElement = elementVector %Busco el los elementos afectados
        
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:);
    
    %Integro por gauss
    for iPoint = 1:nPointInterger
        LocalPointInterger=PointIntegration(iPoint,:);
        
        [Ni,N]=getShapeFunctions(LocalPointInterger,elementType);
        shapeFunctionsDerivatives = getShapeFunctionsDerivatives(LocalPointInterger,elementType);
        ShapeFunction=Ni;
        jacobian = shapeFunctionsDerivatives*elementalNodesPosition;
        
        RadiusNodes = elementalNodesPosition(:,1);
        r =ShapeFunction*RadiusNodes;
        Traction=ones(nElementalNodes,2).*Sigma;     
       
       
        %Integro la carga en los puntos de gauss
        wp=gaussPointsWeight(iPoint);
     
        pointLoadsArray(elementNodesArray(iElement,:),:)=pointLoadsArray(elementNodesArray(iElement,:),:)...
            -ShapeFunction'*ShapeFunction* Traction*r*jacobian(1,1)*wp;
         pointLoadsArray(elementNodesArray(iElement,:),:)=pointLoadsArray(elementNodesArray(iElement,:),:)...
             +ShapeFunction'*ShapeFunction* Traction*r*jacobian(1,2)*wp;
         
    end
    
    %Integro simbolicamente
    Nint = IntegerShapeFunctions(elementType,Side);
    pointLoadsArray2(elementNodesArray(iElement,:),:)=pointLoadsArray2(elementNodesArray(iElement,:),:)...
            +Nint* Traction*r*jacobian(1,1);
      
end