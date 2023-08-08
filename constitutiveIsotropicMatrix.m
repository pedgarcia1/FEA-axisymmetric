function [constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,youngModulus,poissonRatio)
% Constitutive Isotropic Matrix Calculator
%
% constitutiveMatrix:   Constitutive matrix 
%
% problemType: 'Stress' 'Strain' 'Axisymmetric' '3D'
% youngModulus: Young's Modulus
% poissonRatio: Poission's Ratio


switch problemType
    case 'Stress'
        constitutiveMatrix = youngModulus / (1 - poissonRatio^2)*[  1.0             poissonRatio    0.0
                                                                    poissonRatio    1.0             0.0
                                                                    0.0             0.0             (1 - poissonRatio)/2 ];
    case 'Strain'
        constitutiveMatrix = youngModulus / (1 + poissonRatio) /(1 - 2*poissonRatio) * [1.0-poissonRatio    poissonRatio        0.0
                                                                                        poissonRatio        1.0-poissonRatio    0.0
                                                                                        0.0                 0.0                 0.5 - poissonRatio];
    case 'Axisymmetric'
        constitutiveMatrix = youngModulus / (1 + poissonRatio) /(1 - 2*poissonRatio) * [1.0-poissonRatio    poissonRatio         poissonRatio        0.0
                                                                                        poissonRatio        1.0-poissonRatio     poissonRatio        0.0
                                                                                        poissonRatio        poissonRatio         1.0-poissonRatio    0.0
                                                                                        0.0                 0.0                  0.0                 0.5 - poissonRatio];
    case '3D'
        constitutiveMatrix = youngModulus / (1 + poissonRatio) /(1 - 2*poissonRatio) * [1.0-poissonRatio    poissonRatio         poissonRatio        0.0     0.0     0.0
                                                                                        poissonRatio        1.0-poissonRatio     poissonRatio        0.0     0.0     0.0
                                                                                        poissonRatio        poissonRatio         1.0-poissonRatio    0.0     0.0     0.0
                                                                                        0.0                 0.0                  0.0                 0.5 - poissonRatio  0.0                 0.0
                                                                                        0.0                 0.0                  0.0                 0.0                 0.5 - poissonRatio  0.0
                                                                                        0.0                 0.0                  0.0                 0.0                 0.0                 0.5 - poissonRatio];
end