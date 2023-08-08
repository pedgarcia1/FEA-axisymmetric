function [dofVector]=convertNode2Dof(nodesVector,nNodalDof)
% Node numbering to dof converter
% 
% matrixAssembly(elementType,elementNodesArray,nodesPositionArray)
%
% dofVector:        Vector with dof numbering
%
% nodesVector:      Vector with nodal numbering
% nNodalDof:        Number of nodal dof
%

%% Definitions
nNodes=size(nodesVector,2);
dofVector=zeros(1,nNodes*nNodalDof);

% Calculation
for iNodalDof=1:nNodalDof
    dofVector(iNodalDof:nNodalDof:nNodes*nNodalDof)=nNodalDof*nodesVector-(nNodalDof-iNodalDof);
end
