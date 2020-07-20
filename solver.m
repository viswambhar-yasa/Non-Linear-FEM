%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%% This is Solver file 
% Newton Raphson is implemented in this file 
%  Inputs:
%    Parameters  : Input parameters
%    Q           : Modulus parameter for viso-elastic material
%  Outputs:
%       elements        : Contains all the information of the fem analysis like
%                         nodes,displacements,stress
%       Global_disp     : final displacement of the nodes are stored.
%       analyticalsol   : Analytical displacements of each node.
%       strain_history  : contain displacement at each time step of the outer node
%       disp_t          : Contains displacement at diffeent time steps
%       Convergence     : contains number of iteration for solution to converge for each time step
function [element_table,Global_disp,analyticalsol,strain_history,Convergence,disp_t]=solver(parameter,Q)
    dt=1/parameter(11); %time step
    Max_pressure=parameter(7);
    radii=parameter(5);
    nelem=parameter(10);
    tl=parameter(8);
    tf=parameter(9);
    % To implement gauss quadrature the gauss points and weights are declared
    gauss_point=0; 
    gauss_weight=2;
    K_Global=zeros(nelem+1);
    F_internal=zeros(nelem+1,1);
    F_external=zeros(nelem+1,1);
    delta_disp=zeros(nelem+1,1);
    Global_disp=zeros(nelem+1,1);
    nodes=meshGenerator(nelem);
    %Element table:contains all the values for analysis is declared
    element_table=zeros(nelem,12);
    element_table(:,1) = 1:1:nelem;
    element_table(:,2)= nodes(1:end-1);
    element_table(:,3)= nodes(2:end);
    Ke(:,:,nelem)=zeros(2);
    Fe(:,:,nelem)=zeros(2,1);
    disp_t=zeros(nelem+1,2);
    flag=1;
    for T=0:dt:tf
        %Loading scaling is performed  
        if (T<=tl) 
            F_external=zeros(nelem+1,1);
            load_scaling=(1/tl)*T;
            P=Max_pressure*radii*load_scaling;
            F_external(1,1)=P;
        else
            F_external=zeros(nelem+1,1);
            P=Max_pressure*radii;
            F_external(1,1)=P;
        end
        iteration=0;
        %implementation of Newton Raphson method
        while (iteration<20)
            Ke(:,:,nelem)=zeros(2);
            Fe(:,:,nelem)=zeros(2,1);
            K_Global=zeros(nelem+1);
            F_internal=zeros(nelem+1,1);
            F_external=zeros(nelem+1,1);
            delta_disp=zeros(nelem+1,1); %#ok<*PREALL>
            F_external(1,1)=P;
            for i=1:nelem
                %Calculating elemental stiffness and internal forces from elementroutine for each node
                [Ke(:,:,i),Fe(:,:,i),stress,overstress]=elementroutine(element_table(i,:),gauss_point,gauss_weight,T,Q);
                %Assembling Global stiffness and internal force
                K_Global(i:i+1,i:i+1)=K_Global(i:i+1,i:i+1)+Ke(:,:,i);
                F_internal(i:i+1,1)=F_internal(i:i+1,1)+Fe(:,:,i);
                element_table(i,8)=stress(1);
                element_table(i,9)=stress(2);
                element_table(i,10)=overstress(1);
                element_table(i,11)=overstress(2);
            end
            %calculating delta displacements and add it to global displacements
            delta_disp=linsolve(K_Global,(F_external-F_internal));
            Global_disp=Global_disp+delta_disp;
            
            if T==(tl/2) 
                disp_t(:,1)=Global_disp;
            elseif T==tl 
                disp_t(:,2)=Global_disp;    
            elseif T==tf 
                disp_t(:,3)=Global_disp;
            end
            
            
            % the data is stored in element_table to implement NRM in next iteration
            for j=1:nelem
                element_table(j,4)=element_table(j,4)+delta_disp(j);
                element_table(j,5)=element_table(j,5)+delta_disp(j+1);
                element_table(j,6)=delta_disp(j);
                element_table(j,7)=delta_disp(j+1);
            end
 
            % convergence for each loop is calculated and stored, 
            Residual=F_internal-F_external;
            if (max(abs(Residual)))<0.005*(max(abs(F_internal)))|| (max(abs(delta_disp)))<0.005*(max(abs(Global_disp)))
                Convergence(flag,1)=T;
                Convergence(flag,2)=iteration;
                break;
            else  
                iteration=iteration+1;
            end
        end
        %stores strain of the last node for each time step
        strain_history(flag,1)=T;
        strain_history(flag,2)=Global_disp(nelem+1,1);
        strain_history(flag,3)=Global_disp(floor((nelem+1)/2),1);
        strain_history(flag,4)=Global_disp(1,1);
        flag=flag+1;
    end
    analyticalsol=analytical_sol();
end