function [e]=calcErr(phi,h)

phiInf=(phi(1)*h(2)-phi(2)*h(1))/(h(2)-h(1));

e=(phi(2)-phiInf)/phiInf;
% e=phi(2);
% e=phi(2)-phi(1);

end
