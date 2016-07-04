
%% PARAMATERS
ExpData_folder = 'C:\Users\RDickes\Documents\GitHub\PSA_SFERAII_gitHub.git\ExperimentalData\';
ExpData_filename = 'FullDay_2016_06_29';

SimData_folder = 'C:\Users\RDickes\Google Drive\PhD\Modelica\Results\PSA_SFERAII\';
SimData_filename = 'PTTL_SF_basic';


%% EXPERIMENTAL DATA LOADING
load([ExpData_folder ExpData_filename '\' ExpData_filename '.mat'])
ExpData = point;


%% SIMULATION RESULTS LOADING
SimData = dymload([SimData_folder SimData_filename]); 
time_SimData = dymget(SimData,'Time'); 
T_ptc_ex_SimData = dymget(SimData,'SensTex.T'); 
%T_ptc_su_SimData = dymget(SimData,'SensTsu.T'); 
M_dot_SimData = dymget(SimData,'Supply.in_Mdot'); 


%% FIGURE
% time sequence defintion
time_start_plot = 46000;
time_stop_plot = time_SimData(end);

% index evaluation for simulation results
[ ~,index_Sim_start] = min(abs(time_SimData - time_start_plot));
[ ~,index_Sim_stop] = min(abs(time_SimData - time_stop_plot));
vec_Sim_plot = index_Sim_start:index_Sim_stop;

% index evaluation for experimental data
[ ~,index_Exp_start] = min(abs(ExpData.time_sec - time_start_plot));
[ ~,index_Exp_stop] = min(abs(ExpData.time_sec - time_stop_plot));
vec_Exp_plot = index_Exp_start:index_Exp_stop;

% figure creation
LW = 1.5;
figure('units','normalized','outerposition',[0 0 1 1])
hold on
Line (1) = plot(ExpData.time_sec(vec_Exp_plot),ExpData.T_ptc_su(vec_Exp_plot)-273.15, 'LineWidth', LW);
Leg{1} = 'T_{su,ptc}';
Line (2) = plot(time_SimData(vec_Sim_plot),T_ptc_ex_SimData(vec_Sim_plot)-273.15, 'LineWidth', LW);
Leg{2} = 'T_{ex,sim}';
Line (3) = plot(ExpData.time_sec(vec_Exp_plot),ExpData.T_ptc_ex(vec_Exp_plot)-273.15, 'LineWidth', LW);
Leg{3} = 'T_{ex,exp}';
hold off
grid on
ylabel('Outlet temperature [°C]')
xlabel('Day time [sec]')
legend(Line, Leg)





