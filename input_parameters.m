%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%% PARAMETER FILE : This are the input parameter for Visco-Elastic Material
%The time step and number of element shoud be choosen based on the computer performance and convergence
function parameters=input_parameters()
    E=200000;       %Youngs Modulus     1
    v=0.30;         %poisson ratio      2
    Q=100000;       %charatesti modulus 3
    T=1;            %Time char          4
    a=30;           %inner radius       5 
    b=60;           %outer radius       6
    Pressure=140;   %load               7
    ltime=2;        %maxloading time    8
    ftime=10;       %final time         9 

    nelem=30;       %No of elements     10
    noofsteps=10;   % 1/dt              11  

    parameters=[E,v,Q,T,a,b,Pressure,ltime,ftime,nelem,noofsteps];
end