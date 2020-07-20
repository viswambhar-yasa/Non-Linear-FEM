%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%% Material Routine file : To get material stiffness matrix and stress 
%   Inputs:
%       strain          : nodal strains from previous iteration
%       dstrain         : delta strains to calculate overstress
%       pr_overstress   : previous overstress of the nodes
%       T               : current Time step
%       Q               : Modulus to calculate overstress
%   Outputs:
%       Ct          : Material stiffness matrix
%       Stress      : sigma rr and ff for each element
%       overstress  : overstress for current timestep which will be stored in element_table
%REFERENCE : Nonlinear Finite Element Methods by Peter Wriggers
%            3.3.3 Visco-Elastic and Visco-Plastic Material Behaviour
%            6.2.1 Viscoelastic Material

function [Ct,stress,overstress]=materialroutine(strain,dstrain,pr_overstress,T,Q)
    params=input_parameters();
    T=params(4);
    dt=1/params(11);
    E=params(1);
    mu=params(2);
    G=E/((1+mu)*(1-2*mu));
    C=G*[1-mu,mu;mu,1-mu];
    overstress=[0;0]; 
    dev=[0;0];
    %hydrostatic strain
    Hyd_strain = (dstrain(1,1)+dstrain(2,1))/3;
    %deviatoric strain 
    dev(1,1)=dstrain(1,1) -  Hyd_strain;
    dev(2,1)= dstrain(2,1) - Hyd_strain;
    overstress= (1/(1+(dt/(2*T)))) * ((pr_overstress*(1-(dt/(2*T)))) + Q*dev);
    Ct=[0,0;0,0];
    %material stiffness tangent
    Ct(1,1)=C(1,1)+((2/3)*Q/(1+(dt/(2*T))));
    Ct(1,2)=C(1,2)-((1/3)*Q/(1+(dt/(2*T))));
    Ct(2,1)=C(2,1)-((1/3)*Q/(1+(dt/(2*T))));
    Ct(2,2)=C(2,2)+((2/3)*Q/(1+(dt/(2*T))));
    %overall stress is calculated 
    stress= C*strain + overstress;
end