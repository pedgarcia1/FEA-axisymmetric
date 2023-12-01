function [Nint] = IntegerShapeFunctions(elementType,Side)

syms ksi eta;

switch elementType
    
    case 'Q9'
        
        N9 =      (1 - ksi^2)*(1 - eta^2);
        N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
        N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
        N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
        N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
        N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
        N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
        N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
        N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);
        
        N = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
        
        
    case 'AHMAD9'
        
        N9 =      (1 - ksi^2)*(1 - eta^2);
        N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
        N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
        N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
        N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
        N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
        N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
        N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
        N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);
        
        N= [N1 N2 N3 N4 N5 N6 N7 N8 N9];
        
    case 'Q8'
        
        
        N8 = 0.50*(1 - ksi  )*(1 - eta^2);
        N7 = 0.50*(1 - ksi^2)*(1 + eta  );
        N6 = 0.50*(1 + ksi  )*(1 - eta^2);
        N5 = 0.50*(1 - ksi^2)*(1 - eta  );
        N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
        N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
        N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
        N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);
        
        N = [N1 N2 N3 N4 N5 N6 N7 N8];
        
        
    case 'AHMAD8'
        
        N8 = 0.50*(1 - ksi  )*(1 - eta^2);
        N7 = 0.50*(1 - ksi^2)*(1 + eta  );
        N6 = 0.50*(1 + ksi  )*(1 - eta^2);
        N5 = 0.50*(1 - ksi^2)*(1 - eta  );
        N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
        N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
        N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
        N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);
        
        N= [N1 N2 N3 N4 N5 N6 N7 N8];
        
        
    case 'Q4'
        
        N4 = 0.25*(1 - ksi)*(1 + eta);
        N3 = 0.25*(1 + ksi)*(1 + eta);
        N2 = 0.25*(1 + ksi)*(1 - eta);
        N1 = 0.25*(1 - ksi)*(1 - eta);
        
        N = [N1 N2 N3 N4];
        
    case 'AHMAD4'
        
        N4 = 0.25*(1 - ksi)*(1 + eta);
        N3 = 0.25*(1 + ksi)*(1 + eta);
        N2 = 0.25*(1 + ksi)*(1 - eta);
        N1 = 0.25*(1 - ksi)*(1 - eta);
        
        N= [N1 N2 N3 N4];
    otherwise
        disp('Unknown method.')
end

switch Side
    case 'Up'
        Ninteger=int(transpose(N)*N,ksi,[-1 1]);
        eta=1;
        Nint=eval(Ninteger);
        
    case 'Down'
        Ninteger=int(transpose(N)*N,ksi,[-1 1]);
        eta=-1;
        Nint=eval(Ninteger);
    case 'Left'
        Ninteger=int(transpose(N)*N,eta,[-1 1]);
        ksi=-1;
        Nint=eval(Ninteger);
    case 'Right'
        Ninteger=int(transpose(N)*N,eta,[-1 1]);
        ksi=1;
        Nint=eval(Ninteger);
    otherwise
        disp('Unknown method.')
end



end


