%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%% This element routine file
% performs calculation to get nodal stiffness and internal force
% Inputs:   
%       element         : contains all info of the nodes from previous iteration
%       gauss_points    : [0  0]
%       gauss_weight    : 0
%       T               : Time increment required in material routine
%       Q               : Modulus required in material routine
% Outputs:
%       Kt    : nodal stiffness matrix
%       Fe    : nodal internal force
%       Stress: nodal stress rr and ff

function [Kt,Fe,stress,overstress]=elementroutine(element,gauss_point,gauss_weight,T,Q)
        %nodal position ,du and displacements for previous iteration are extracted
        R=[element(2),element(3)];
        U=[element(4);element(5)];
        DU=[element(6);element(7)];
        pr_overstress=[element(10);element(11)];
        %Jacobian and Jacobian_inverse are calculated
        J=(R(2)-R(1))/2;
        inve_J=1/J;
        %shape function and strain-displacement are calculated
        N=shape_function(gauss_point);
        B=B_Mat(inve_J,R,gauss_point);
        B_Trans=B';
        %Inputs for Material routine are calculated
        strain=B*U;
        dstrain=B*DU;
        NR=N(1)*R(1)+N(2)*R(2);
        %material stiffness matrix is obtained from material routine
        [C,stress,overstress]=materialroutine(strain,dstrain,pr_overstress,T,Q);
        % Nodal stiffness and internal force are calculated
        Fe=gauss_weight*B_Trans*stress*NR*J;
        Kt=gauss_weight*B_Trans*C*B*NR*J;     
end
% Shape function calculated for each gauss points
% Inputs:
%   Gauss_points :[0 0]
function N=shape_function(gauss_point)
    N1=0.5*(1-gauss_point);
    N2=0.5*(1+gauss_point);
    N=[N1,N2];
end

% Strain-displacement matrix 
%    Inputs:
%           inve_J  : Jacobain Inverse
%           R       : nodal vector [r1 r2]
%           Eplison :gauss points [0 0]
function B_Mat=B_Mat(inve_J,R,eplison)
    N=shape_function(eplison);
    B_Mat(1,1)=-0.5*inve_J;
    B_Mat(1,2)=0.5*inve_J;
    B_Mat(2,1)=N(1)/(N(1)*R(1)+N(2)*R(2));
    B_Mat(2,2)=N(2)/(N(1)*R(1)+N(2)*R(2));
end