%% Author               : Yasa Viswambhar Reddy
%% Matriculation number : 65074
%    Non-Linear Finite element method for axisymmertic, plane strain for a viscoelastic material
%    This is the main file which get inputs and calls solver and plots results
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2020
%    lecturer in charge: Dr. Geralf Hütter
clc;
clear;
tic;% to calculate time taken for the process to run 
parameter=input_parameters();
nodes=meshGenerator(parameter(10));
prompt = '[Y]-Non linear [N]-linear (Q=0) :';
str = input(prompt,'s');
if str=='N'
    Q=0;
    [element_Table,Global_disp,analyticalsol,strain_history,Convergence,disp_t]=solver(parameter,Q);
    %plotting the results
    figure(1);
        plot(nodes,analyticalsol,'r-.o',nodes,Global_disp,'b:*');
        title('Analytical solution vs Numerical solution at final time')
        xlabel('Radius r (mm)');
        ylabel('Radial displacement u_{r} (mm)');
        legend('Analytical','Numerical');
        legend('boxoff');
        grid on;
        grid minor;
        print('Linear','-dpng')
    figure(2);
        plot(Convergence(:,1),Convergence(:,2));
        xlim([0 10]);
        ylim([0 2]);
        title('Convergence at each time step');
        xlabel('Time (s) ');
        ylabel('Number of iterations');
        grid on;
        grid minor;
        print('LinearConvergence','-dpng');
else
    Q=parameter(3);
    [element_Table,Global_disp,analyticalsol,strain_history,Convergence,disp_t]=solver(parameter,Q);
    stress_rr=element_Table(:,8);
    stress_fi=element_Table(:,9);
    %plotting the results
    figure(1);
        plot(nodes,analyticalsol,'r-.o',nodes,Global_disp,'b:*');
        title('Analytical solution vs Numerical solution')
        xlabel('Radius r (mm)');
        ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
        legend('Analytical','Numerical');
        legend('boxoff');
        grid on;
        grid minor;
    print('AnalyticalNumerical','-dpng')
    figure(2);
        %plot(strain_history(:,1),strain_history(:,2),strain_history(:,1),strain_history(:,3),strain_history(:,1),strain_history(:,4),'linewidth',1.5);
        plot(strain_history(:,1),strain_history(:,2),strain_history(:,1),strain_history(:,3),'linewidth',1.5);
        title('Time history of the widening of the pipe u_{r} at b and a','Interpreter', 'tex')
        xlabel('Time t (s)');
        ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
        legend('last node b ','first node a','location','SouthEast');
        legend('boxoff'); 
        grid on;
        grid minor;
    print('wideningpipe','-dpng')
    figure(3);
        subplot(2,1,1);
            plot(nodes(1:end-1),stress_rr,'b-*','linewidth',1);
            set(gcf,'Position',[510 100 500 500])
            title('{\sigma}_{r r} vs  r','Interpreter', 'tex');
            xlabel('radius r (mm)');
            ylabel('{\sigma}_{r r}(MPa)','Interpreter', 'tex');
            grid on;
            grid minor;
        subplot(2,1,2);
            plot(nodes(1:end-1),stress_fi,'r-o','linewidth',1);
            set(gcf,'Position',[510 100 500 500])
            title('{\sigma}_{\phi \phi} vs  r','Interpreter', 'tex');
            xlabel('radius r (mm)');
            ylabel('{\sigma}_{\phi \phi} (MPa)','Interpreter', 'tex');
            grid on;
            grid minor;
    print('stressrrff','-dpng')
    figure(4);
        plot(Convergence(:,1),Convergence(:,2));
        xlim([0 10]);
        ylim([0 20]);
        title('Convergence at each time step');
        xlabel('Time (s) ');
        ylabel('Number of iterations');
        grid on;
        grid minor;
    print('NonLinearConvergence','-dpng');
    figure(5);
        plot(nodes,disp_t(:,1),'g--o',nodes,disp_t(:,2),'r--o',nodes,disp_t(:,3),'b:o',nodes,Global_disp,'k-.*');
        title('displacements at different time steps')
        xlabel('Radius r (mm)');
        ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
        legend('Half time to max load time','Max load time','final time','Analytical Solution' );
        legend('boxoff');
        grid on;
        grid minor;
    print('displacementtNL','-dpng')
end
%Generates a csv file from the matrix 
header = {'ELEMENTS','Radius (a)','Radius(b)','disp u(a)','disp u(b)','delta_disp du(a)','delta_disp du(b)','stress rr','stress ff','overstress rr','overstress ff','null'};
ELEMENT_TABLES = [header; num2cell(element_Table)];
writecell(ELEMENT_TABLES,'ELEMENT_TABLE.csv')
toc;