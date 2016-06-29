
% Code initializer
clear all
close all
clc

%% PARAMETERS
filename = '27-06-2016_bis.txt';
folder_path = [cd '\ExperimentalData\' ];

%% DATA IMPORTATION
formatSpec = '%{dd/MM/yyyy}D%{HH:mm:ss}D%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';	
Data = readtable([ folder_path filename],'Delimiter','\t', 'Format',formatSpec);


%% DATA POSTREATEMENT
time_global_start = Data.Hora(1);
time_global_stop = Data.Hora(end);

% Sampling time - definition
time_sample_start = time_global_start; %'01:00:00';
time_sample_stop = time_global_stop; %'02:00:00';

[ ~,index_sample_start] = min(abs(Data.Hora - time_sample_start));
[ ~,index_sample_stop] = min(abs(Data.Hora - time_sample_stop));
vec_sample = index_sample_start:index_sample_stop;



%% FIGURES
LW = 1;
LS = '-';

Vector2Plot =   [1,         2,          3,          4,          5,          6,          7,          8,          9,          10,         11,         12,         13,         14,                   15,               16];
Variable2Plot = {'TA013',   'TA015',    'TA016',    'TA030',    'TA031',    'TA060',    'TA066',    'TA011',    'TA012',    'FA032',    'FA023',    'PA020',    'PA052',    'pocision_EURO',      'consigna_EURO',  'Incidencia'};
Label2Plot =    {'TA013',   'TA015',    'TA016',    'TA030',    'TA031',    'TA060',    'TA066',    'TA011',    'TA012',    'FA032',    'FA023',    'PA020',    'PA052',    'pocision_EURO',      'consigna_EURO',  'Incidencia'};

% Global results - TEMPERATURE PROFILES
figure('Name', 'Temperature Profile - Global')
hold all
j= 0;
for k_Tplot = 1:9
    j = j+1;
    eval(['Tplot(j) = plot(Data.Hora, Data.' Variable2Plot{k_Tplot} ', ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegT{j} = Label2Plot{k_Tplot};
end
hold off
legend(Tplot, LegT)
grid on
xlabel('Time')
ylabel('Temperature [°C]')

% Global results - FLOW PROFILES
figure('Name', 'Flow Rate Profile - Global')
hold all
j= 0;
for k_Fplot = 10:11
    j = j+1;
    eval(['Fplot(j) = plot(Data.Hora, Data.' Variable2Plot{k_Fplot} ', ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegF{j} = Label2Plot{k_Fplot};
end
hold off
legend(Fplot, LegF)
grid on
xlabel('Time')
ylabel('Flow rate [kg/s - m3/s]')

% Global results - PRESSURE PROFILES
figure('Name', 'Pressure Profile - Global')
hold all
j= 0;
for k_Pplot = 12:13
    j = j+1;
    eval(['Pplot(j) = plot(Data.Hora, Data.' Variable2Plot{k_Pplot} ', ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegP{j} = Label2Plot{k_Pplot};
end
hold off
legend(Pplot, LegP)
grid on
xlabel('Time')
ylabel('Temperature [°C]')


% Sample results - TEMPERATURE PROFILES
figure('Name', 'Temperature Profile - Sample')
hold all
j= 0;
for k_Tplot = 1:9
    j = j+1;
    eval(['Tplot(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k_Tplot} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegT{j} = Label2Plot{k_Tplot};
end
hold off
legend(Tplot, LegT)
grid on
xlabel('Time')
ylabel('Temperature [°C]')

% Sample results - FLOW PROFILES
figure('Name', 'Flow Rate Profile - Sample')
hold all
j= 0;
for k_Fplot = 10:11
    j = j+1;
    eval(['Fplot(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k_Fplot} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegF{j} = Label2Plot{k_Fplot};
end
hold off
legend(Fplot, LegF)
grid on
xlabel('Time')
ylabel('Flow rate [kg/s - m3/s]')

% Sample results - PRESSURE PROFILES
figure('Name', 'Pressure Profile - Sample')
hold all
j= 0;
for k_Pplot = 12:13
    j = j+1;
    eval(['Pplot(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k_Pplot} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
    LegP{j} = Label2Plot{k_Pplot};
end
hold off
legend(Pplot, LegP)
grid on
xlabel('Time')
ylabel('Temperature [°C]')

%% EXPORT RESULTS : 
if 0
    clear point
    Point_name = 'EXP1_27_06_2016';
    vec_export = vec_sample;
    point.vector_sample = vec_sample;
    point.T_ptc_su = Data.TA060(vec_sample); %or TA058 or TA046?
    point.T_ptc_ex = Data.TA066(vec_sample); % or TA074 or TA076?
    point.P_ptc_su = Data.PA052(vec_sample); % bar
    point.P_tank = Data.PA021(vec_sample);% bar
    point.P_pp_ex = point.PA020; %bar
    point.M_dot_htf = Data.FA032(vec_sample);% kg/s
    point.V_dot_htf = Data.FA023(vec_sample);% m3/s
    point.DNI  = Data.IA028(vec_sample);% W/m2
    point.T_amb = Data.TA029(vec_sample);% 
    point.theta = Data.Hora(vec_sample);% or???
    point.V_wind_5 = Data.ST082(vec_sample);% 
    point.V_wind_12 = Data.ST072(vec_sample);% 
    point.D_wind_5 = Data.WD089(vec_sample);% compared to the north (+90 is east)
    point.D_wind_12 = Data.WD088(vec_sample);% compared to the north (+90 is east)
    point.time_day = Data.Hora(vec_sample);
    point.time_sec = NaN*ones(length(vec_sample),1);
    point.time_vs_DNI_text = '';
    for k = 1:length(vec_sample)
        [Y,M,D,H,MN,S] = datevec(Data.Hora(vec_sample(k))-Data.Hora(vec_sample(1)));
        point.time_sec(k,1) = (((Y*365 + M*30 + D)*24 + H)*60 + MN)*60 + S;        
        point.time_vs_DNI_text = [point.time_vs_DNI_text num2str(point.time_sec(k,1)) ',' num2str(point.DNI(k,1)) '; '];
        point.time_vs_Tamb_text = [point.time_vs_Tamb_text num2str(point.time_sec(k,1)) ',' num2str(point.T_amb(k,1)) '; '];
        point.time_vs_theta_text = [point.time_vs_theta_text num2str(point.time_sec(k,1)) ',' num2str(point.theta(k,1)) '; '];
        point.time_vs_Mdot_text = [point.time_vs_Mdot_text num2str(point.time_sec(k,1)) ',' num2str(point.M_dot_htf(k,1)) '; '];
        point.time_vs_Tsu_text = [point.time_vs_Tsu_text num2str(point.time_sec(k,1)) ',' num2str(point.T_ptc_su(k,1)) '; '];
        point.time_vs_Psu_text = [point.time_vs_Psu_text num2str(point.time_sec(k,1)) ',' num2str(point.P_ptc_su(k,1)) '; '];
        point.time_vs_Pex_text = [point.time_vs_Pex_text num2str(point.time_sec(k,1)) ',' num2str(point.P_ptc_su(k,1)) '; '];

    end
    
    eval([ Point_name ' = point;'])
    eval(['save(''' Point_name '.mat'', ''' Point_name ''' ) '])
end