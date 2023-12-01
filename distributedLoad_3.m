function [pointLoadsArray] = distributedLoad_3(elementType,Lside,direc,pointLoadsArray,nodesPositionArray,pressureNormal)

q = @(r) pressureNormal*2*pi*r;

switch elementType
    case 'Q8'

        matQ = @(a) a/15*[4 2 -1;
            2 16 2;
            -1 2 4];

        for ni = 1:2:length(Lside)-1

            a(1) = norm(nodesPositionArray(Lside(ni),:) - nodesPositionArray(Lside(ni+1),:));
            a(2) = norm(nodesPositionArray(Lside(ni+1),:) - nodesPositionArray(Lside(ni+2),:));
            a = mean(a);

            q_vec(1) = q(nodesPositionArray(Lside(ni),1));
            q_vec(2) = q(nodesPositionArray(Lside(ni+1),1));
            q_vec(3) = q(nodesPositionArray(Lside(ni+2),1));

            Fi = matQ(a)*q_vec';

            pointLoadsArray(Lside(ni:(ni+2)),direc) = Fi + pointLoadsArray(Lside(ni:(ni+2)),direc);

        end

    case 'Q4'

        matQ = @(a) a/3*[2 1;1 2];

        for ni = 1:length(Lside)-1

            a = 1/2*norm(nodesPositionArray(Lside(ni),:) - nodesPositionArray(Lside(ni+1),:));

            q_vec(1) = q(nodesPositionArray(Lside(ni),1));
            q_vec(2) = q(nodesPositionArray(Lside(ni+1),1));

            Fi = matQ(a)*q_vec';

            pointLoadsArray(Lside(ni:(ni+1)),direc) = Fi + pointLoadsArray(Lside(ni:(ni+1)),direc);

        end
    case 'Q4Esfera'

    case 'Q8int'

        

end
end

