function [PointIntegration,gaussPointsWeight] =getPointInteger(Side,elementType)
switch elementType
    case 'Q4'
        nGaussPoints = 2;
        elementShape='Line';
        [gaussPointsLocation,gaussPointsWeight] = getGaussPoints(elementShape,nGaussPoints);
        switch Side
            case 'Up'
                  PointIntegration=[gaussPointsLocation',[1;1]];
            case 'Down'
                 PointIntegration=[gaussPointsLocation',[-1;-1]];
            case 'Left'
                 PointIntegration=[[-1;-1],gaussPointsLocation'];
            case 'Right'
                 PointIntegration=[[1;1],gaussPointsLocation'];
        end
    case 'Q8'
        nGaussPoints = 3;
        elementShape='Line';
        [gaussPointsLocation,gaussPointsWeight] = getGaussPoints(elementShape,nGaussPoints);
        switch Side
            case 'Up'
                PointIntegration=[gaussPointsLocation',[1;1;1]];
            case 'Down'
               PointIntegration=[gaussPointsLocation',[-1;-1;-1]];
            case 'Left'
                PointIntegration=[[-1;-1;-1],gaussPointsLocation'];
            case 'Right'
               PointIntegration=[[1;1;1],gaussPointsLocation'];
        end
        
end
