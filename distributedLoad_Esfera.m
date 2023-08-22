function [pointLoadsArray] = distributedLoad_Esfera(elementType,Lside,pointLoadsArray,nodesPositionArray,pressureNormal,centro)

q = @(r) pressureNormal*2*pi*r;

switch elementType
    case 'Q8'

    case 'Q4'

        matQ = @(a) a/3*[2 1;1 2];

        

        for ni = 1:length(Lside)-1

            a = 1/2*norm(nodesPositionArray(Lside(ni),:) - nodesPositionArray(Lside(ni+1),:));
            
            ang(1) = atan2d(abs(centro(2)-nodesPositionArray(Lside(ni),2)),abs(centro(1)-nodesPositionArray(Lside(ni),1)));
            ang(2) = atan2d(abs(centro(2)-nodesPositionArray(Lside(ni+1),2)),abs(centro(1)-nodesPositionArray(Lside(ni+1),1)));

            % cargas en R
            q_vec(1) = cosd(ang(1))*q(nodesPositionArray(Lside(ni),1));
            q_vec(2) = cosd(ang(2))*q(nodesPositionArray(Lside(ni+1),1));

            Fi = matQ(a)*q_vec';

            pointLoadsArray(Lside(ni:(ni+1)),1) = Fi + pointLoadsArray(Lside(ni:(ni+1)),1);

            % cargas en Z
            q_vec(1) = sind(ang(1))*q(nodesPositionArray(Lside(ni),1));
            q_vec(2) = sind(ang(2))*q(nodesPositionArray(Lside(ni+1),1));

            Fi = matQ(a)*q_vec';

            pointLoadsArray(Lside(ni:(ni+1)),2) = Fi + pointLoadsArray(Lside(ni:(ni+1)),2);

        end
    case 'Q4Esfera'

end
end

