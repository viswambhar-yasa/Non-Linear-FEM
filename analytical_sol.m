%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%% Anaytical Solution is calculated using this function for each node
function analyticalsol=analytical_sol()
    parameter=input_parameters();
    nelem=parameter(10);%Number of elements
    Rin=parameter(5); %Inner Radius
    Rout=parameter(6); %outer Radius
    E=parameter(1); %Youngs Modulus
    v=parameter(2); %poissions ratio
    maxPres=parameter(7); %Pressure
    Constant=(1+v)*(maxPres/E)*((Rin*Rin)/((Rout*Rout)-(Rin*Rin)));
    r=meshGenerator(nelem);
    analyticalsol=zeros(nelem+1,1);
    for i=1:nelem+1
        analyticalsol(i)=Constant*(((1-2*v)*r(i))+((Rout*Rout)/r(i)));
    end
end