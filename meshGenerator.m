%% Generate list of position of nodes according to a geometric series
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2020
%    lecturer in charge: Dr. Geralf Hütter
function rnodes=meshGenerator(nelem)
    parameters=input_parameters();
    b=parameters(6); %outer radius
    a=parameters(5); %inner radius
%number of elements
    meshrefinementfactor=2; %ratio of element sizes at outer and inner radius

%ratio between element sizes of subsequent elements for a geometric series
    q=meshrefinementfactor^(1./(nelem-1));
%size of first interval
    dr=(b-a)*(1-q)/(1-meshrefinementfactor*q);
    rnode=a;
    rnodes=[a];
%loop over all elements
for i=1:nelem
        rnode=rnode+dr;
        rnodes=[rnodes;rnode];
        dr=dr*q;
end
end