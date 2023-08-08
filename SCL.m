function [Svm_m,Svm_b,Svm_F,Sfem,Sm,Sb,SF,Sfemx]=SCL(elementStressAtNodes,msh,SCLnod)

SCLstress=zeros(length(SCLnod),3);
% msh.elem.nod(:,1);


for i=1:length(SCLnod)
  [r,c]=find(msh.elem.nod==SCLnod(i));
  SCLelem{i,1}=r;
  SCLelem{i,2}=c;
end
for i=1:length(SCLnod)
SCLstress(i,1)=mean(diag(elementStressAtNodes(SCLelem{i,1},SCLelem{i,2},1)));
SCLstress(i,2)=mean(diag(elementStressAtNodes(SCLelem{i,1},SCLelem{i,2},2)));
SCLstress(i,3)=mean(diag(elementStressAtNodes(SCLelem{i,1},SCLelem{i,2},3)));

end
% dnod=vecnorm((msh.cord(SCLnod,:))')';
dnod=msh.cord(SCLnod,2);
ddnod=dnod-dnod(1);
% SCLstress(1,:)./elementStressAtNodes(6,2)
Sm=trapz(ddnod,SCLstress)/ddnod(end);
Sb=trapz(ddnod,SCLstress.*(ddnod(end)/2-ddnod))*6/ddnod(end)^2;
SF=SCLstress(end,:)-(Sm-Sb);
Svm_m=sqrt(Sm(1)^2+Sm(2)^2+3*Sm(3)^2);
Svm_b=sqrt(Sb(1)^2+Sb(2)^2+3*Sb(3)^2);
Sfem=sqrt(SCLstress(:,1).^2+SCLstress(:,2).^2+3*SCLstress(:,3).^2);
Svm_F=Sfem(end)-(Svm_m+Svm_b);%sqrt(SF(1)^2+SF(2)^2+3*SF(3)^2);
Sfemx=SCLstress(:,1);
% Sm=Sm';Sb=Sb';SF=SF';
end

