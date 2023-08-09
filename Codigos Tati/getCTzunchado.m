function [nCT,CTMatrix] = getCTzunchado(nodeDofs,constraintsRelations,precondCT,nDofTot_U)

nCT=size(constraintsRelations,1);
CTMatrix=zeros(nCT,nDofTot_U);


for iCT=1:nCT
    CTMatrix(iCT,[nodeDofs(constraintsRelations(iCT,1),1) nodeDofs(constraintsRelations(iCT,2),1)])=[-1 1]*precondCT; 
end

end


% function [ CTFrac ,nCREqFrac ] = getCTFrac(nodeDofs,constraintsRelations,precondCT,nDofTot_U)
% %% - Borde de fractura sin interseccion - 
% nCREqFracSinInt    = sum(sum(constraintsRelations>0,2)<4)*3;
% nCRFracSinInt      = sum(sum(constraintsRelations>0,2)<4);
% CTFracSinInt = zeros(nCREqFracSinInt,nDofTot_U);
% constraintsRelationsSinInt = constraintsRelations(sum(constraintsRelations>0,2)<4,:);
% iFila = 1;
% 
% for iCR = 1:nCRFracSinInt
%     midNode = constraintsRelationsSinInt(iCR,1);
%     slave1 = constraintsRelationsSinInt(iCR,2);
% 
% 
%     % EQ FOR v (direccion Y) DOFS %
%     CTFracSinInt(iFila,[nodeDofs(midNode,2) nodeDofs(slave1,2)]) = precondCT*[-1 1];
%     iFila = iFila + 1;
% 
% 
% end
