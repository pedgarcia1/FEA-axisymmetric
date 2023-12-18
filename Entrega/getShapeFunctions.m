function [Ni,N] = getShapeFunctions(evaluationPointsArray,elementType)


nPoints = size(evaluationPointsArray,1);

switch elementType

    case 'Q9'
        N  = zeros(2,18,nPoints);
        Ni = zeros(1, 9,nPoints);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N9 =      (1 - ksi^2)*(1 - eta^2);
            N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
            N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
            N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
            N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);

            Ni(1,:     ,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (1,1:2:17,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (2,2:2:18,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
        end

    case 'AHMAD9'
        N  = zeros(3,3*9,nPoints);
        Ni = zeros(1,  9,nPoints);
        Id = eye(3);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N9 =      (1 - ksi^2)*(1 - eta^2);
            N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
            N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
            N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
            N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);

            Ni(1,:,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (:,:,iPoints) = [ N1*Id, N2*Id, N3*Id, ...
                N4*Id, N5*Id, N6*Id, ...
                N7*Id, N8*Id, N9*Id ];
        end

    case 'Q8'
        N  = zeros(2,16,nPoints);
        Ni = zeros(1, 8,nPoints);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N8 = 0.50*(1 - ksi  )*(1 - eta^2);
            N7 = 0.50*(1 - ksi^2)*(1 + eta  );
            N6 = 0.50*(1 + ksi  )*(1 - eta^2);
            N5 = 0.50*(1 - ksi^2)*(1 - eta  );
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);

            Ni(1,:     ,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (1,1:2:15,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (2,2:2:16,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8];
        end

    case 'AHMAD8'
        N  = zeros(3,3*8,nPoints);
        Ni = zeros(1,  8,nPoints);
        Id = eye(3);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N8 = 0.50*(1 - ksi  )*(1 - eta^2);
            N7 = 0.50*(1 - ksi^2)*(1 + eta  );
            N6 = 0.50*(1 + ksi  )*(1 - eta^2);
            N5 = 0.50*(1 - ksi^2)*(1 - eta  );
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);

            Ni(1,:,iPoints) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (:,:,iPoints) = [ N1*Id, N2*Id, N3*Id, ...
                N4*Id, N5*Id, N6*Id, ...
                N7*Id, N8*Id ];
        end

    case 'Q4'
        N  = zeros(2,8,nPoints);
        Ni = zeros(1,4,nPoints);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N4 = 0.25*(1 - ksi)*(1 + eta);
            N3 = 0.25*(1 + ksi)*(1 + eta);
            N2 = 0.25*(1 + ksi)*(1 - eta);
            N1 = 0.25*(1 - ksi)*(1 - eta);

            Ni(1,:    ,iPoints) = [N1 N2 N3 N4];
            N (1,1:2:7,iPoints) = [N1 N2 N3 N4];
            N (2,2:2:8,iPoints) = [N1 N2 N3 N4];
        end

    case 'AHMAD4'
        N  = zeros(3,3*4,nPoints);
        Ni = zeros(1,  4,nPoints);
        Id = eye(3);
        for iPoints = 1:nPoints
            ksi = evaluationPointsArray(iPoints,1);
            eta = evaluationPointsArray(iPoints,2);

            N4 = 0.25*(1 - ksi)*(1 + eta);
            N3 = 0.25*(1 + ksi)*(1 + eta);
            N2 = 0.25*(1 + ksi)*(1 - eta);
            N1 = 0.25*(1 - ksi)*(1 - eta);

            Ni(1,:,iPoints) = [N1 N2 N3 N4];
            N (:,:,iPoints) = [ N1*Id, N2*Id, N3*Id, N4*Id ];
        end

end

