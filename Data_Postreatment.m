
% Code initializer
clear all
close all
clc

%% PARAMETERS
filename = '29-06-2016_bis.txt';
folder_path = [cd '\ExperimentalData\' ];

%% DATA IMPORTATION
formatSpec = '%{dd/MM/yyyy}D%{HH:mm:ss}D%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
Data = readtable([ folder_path filename],'Delimiter','\t', 'Format',formatSpec);

%% SENSOR NOMENCLATURE
%Fecha      --> Date
%Hora       --> Local day time
%Incidencia	--> theta angle
%TA030      --> T_su_heater
%TA031      --> T_ex_heater
%FA032      --> M_dot measured by Coriolis flow meter
%TA007      --> T_tank_low
%TA008      --> T_tank_med
%TA009      --> T_tank_high
%TA011      --> T_su_cooler
%TA012      --> T_ex_cooler
%TA013      --> T_ex_pump
%TA015      --> T_su_heater2
%TA016      --> T_ex_heater2
%TA017      --> Temperature of the pump
%TA019      --> USELESS
%PA020      --> Pressure supply of the pump
%PA021      --> Pressure in the tank
%CA022      --> Power in the heater
%FA023      --> Volumetric flow rate of the pump
%VA025      --> Load of the pump
%SA026      --> USELESS (old  wind speed)
%DA027      --> USELESS (old wind direction)
%IA028      --> DNI
%TA029	    --> Ambient temperature
%posicion_ALBIASA	---> USELESS
%consigna_ALBIASA	--- USELESS
%posicion_EURO	--> Actual position of the Eurotrough PTC
%consigna_EURO	--> Required position of the Eurotrough PTC
%posicion_URSSA	--> USELESS
%consigna_URSSA	--> USELESS
%TA046      --> USELESS (T_su_eurotrough, but not used)
%TA048      --> USELESS (T_su_albiasa, but not used)
%EMPTY      --> USELESS
%EMPTY      --> USELESS
%TA050      --> USELESS
%EMPTY      --> USELESS
%EMPTY      --> USELESS
%PA052      --> Pressure entrance Eurotrough collectors
%PA053      --> Pressure difference on Eurotrough collectors
%TA054      --> USELESS 
%FA056      --> USELESS
%TA058      --> USELESS
%TA060      --> Supply temperature Eurotrough
%TA066      --> Exhaust temperature Eurotrough
%ST087      --> Wind speed at 5m
%WD088      --> Wind direction at 12m
%ST072      --> Wind speed at 12m
%TA074      --> USELESS
%TA076      --> USELESS
%TA078      --> USELESS
%TA080      --> USELESS
%TA082      --> USELESS
%WD089      --> Wind direction at 5m
%EMPTY      --> USELESS
%EMPTY      --> USELESS
%FA065      --> USELESS
%TA062      --> USELESS
%TA069      --> USELESS
%TA086      --> USELESS
%TA068      --> USELESS
%FA090      --> Flow after the heater
%FA106      --> Flow in Eurotrough (calculated)

%% DATA POSTREATEMENT
time_global_start = '08:00:00';
time_global_stop = Data.Hora(end);
[ ~,index_global_start] = min(abs(Data.Hora - time_global_start));
[ ~,index_global_stop] = min(abs(Data.Hora - time_global_stop));
vec_global = index_global_start:index_global_stop;

% Sampling time - definition
time_sample_start = time_global_start; %'11:30:00';
time_sample_stop = time_global_stop; %'12:55:00';
[ ~,index_sample_start] = min(abs(Data.Hora - time_sample_start));
[ ~,index_sample_stop] = min(abs(Data.Hora - time_sample_stop));
vec_sample = index_sample_start:index_sample_stop;



%% FIGURES
LW = 1;
LS = '-';
plot_all = 1;
plot_sample = 1;

Vector2Plot =   [1,         2,              3,              4,          5,          6,         	7,              8,                	9,                  10,             11,         12,         13,         14,         15,               16              17          18];
Variable2Plot = {'TA029',	'TA060',    	'TA066',        'FA032',	'FA023',    'PA021',    'PA052',        'pocision_EURO', 	'consigna_EURO',    'Incidencia'    'IA028'     'ST087'     'ST072'     'WD089'     'WD088'};
Label2Plot =    {'T_{amb}', 'T_{ptc,su}',	'T_{ptc,ex}',   'M_{dot}',	'V_{dot}',	'P_{tk}',	'P_{ptc,su}',	'pocision_EURO',  	'consigna_EURO',    'Incidencia'    'DNI'       'V_{wd,5}'  'V_{wd,12}'	'D_{wd,5}'  'D_{wd,12}'};


if plot_all
    % Global results - TEMPERATURE PROFILES
    figure('Name', 'Temperature Profile - Global')
    subplot(2,3,1)
    hold all
    j= 0;
    for k = 2:3
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Line Leg k 
    
    % Global results - FLOW PROFILES
    subplot(2,3,2)
    hold all
    j= 0;
    for k = 4:5
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Line Leg k
    
    % Global results - PRESSURE PROFILES
    subplot(2,3,3)
    hold all
    j= 0;
    for k = 6:7
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time') 
    clear Line Leg k
    
    subplot(2,3,4)
    hold all
    j= 0;
    for k = 11
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Leg Line k
    
    subplot(2,3,5)
    hold all
    j= 0;
    for k = 10
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Leg Line k
    
    subplot(2,3,6)
    hold all
    j= 0;
    for k = 12
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_global), Data.' Variable2Plot{k} '(vec_global), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    
    tightfig
end

if plot_sample
    figure('Name', 'Profile - Sample')
    subplot(2,3,1)
    hold all
    j= 0;
    for k = 2:3
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Line Leg k 
    
    % Global results - FLOW PROFILES
    subplot(2,3,2)
    hold all
    j= 0;
    for k = 4:5
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Line Leg k
    
    % Global results - PRESSURE PROFILES
    subplot(2,3,3)
    hold all
    j= 0;
    for k = 6:7
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time') 
    clear Line Leg k
    
    subplot(2,3,4)
    hold all
    j= 0;
    for k = 11
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Leg Line k
    
    subplot(2,3,5)
    hold all
    j= 0;
    for k = 10
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    clear Leg Line k
    
    subplot(2,3,6)
    hold all
    j= 0;
    for k = 12
        j = j+1;
        eval(['Line(j) = plot(Data.Hora(vec_sample), Data.' Variable2Plot{k} '(vec_sample), ''LineStyle'', LS, ''LineWidth'', LW);'])
        Leg{j} = Label2Plot{k};
    end
    hold off
    legend(Line, Leg)
    grid on
    xlabel('Time')
    
    tightfig
end
%% EXPORT RESULTS :
if 1
    clear point
    Point_name = 'FullDay_29_06_2016';
    vec_export = vec_sample;
    point.vector_sample = vec_sample;
    point.T_ptc_su = Data.TA060(vec_sample); 
    point.T_ptc_ex = Data.TA066(vec_sample); 
    point.P_ptc_su = Data.PA052(vec_sample)*1e5; % bar
    point.M_dot_htf = Data.FA032(vec_sample);% kg/s
    point.DNI  = Data.IA028(vec_sample);% W/m2
    point.T_amb = Data.TA029(vec_sample);%
    point.theta = Data.Incidencia(vec_sample)*pi/180;% 
    point.V_wind_5 = Data.ST087(vec_sample);%
    point.D_wind_5 = Data.WD089(vec_sample);% compared to the north (+90 is east)
    point.time_day = Data.Hora(vec_sample);
    point.time_sec = NaN*ones(length(vec_sample),1);
    point.time_vs_DNI_text = '';
    point.time_vs_Tamb_text = '';
    point.time_vs_theta_text = '';
    point.time_vs_Mdot_text = '';
    point.time_vs_Tsu_text = '';
    point.time_vs_Psu_text = '';
    for k = 1:length(vec_sample)
        [Y,M,D,H,MN,S] = datevec(Data.Hora(vec_sample(k))-Data.Hora(vec_sample(1)));
        point.time_sec(k,1) = (((Y*365 + M*30 + D)*24 + H)*60 + MN)*60 + S;
        point.time_vs_DNI_text = [point.time_vs_DNI_text num2str(point.time_sec(k,1)) ',' num2str(point.DNI(k,1)) '; '];
        point.time_vs_Tamb_text = [point.time_vs_Tamb_text num2str(point.time_sec(k,1)) ',' num2str(point.T_amb(k,1)) '; '];
        point.time_vs_theta_text = [point.time_vs_theta_text num2str(point.time_sec(k,1)) ',' num2str(point.theta(k,1)) '; '];
        point.time_vs_Mdot_text = [point.time_vs_Mdot_text num2str(point.time_sec(k,1)) ',' num2str(point.M_dot_htf(k,1)) '; '];
        point.time_vs_Tsu_text = [point.time_vs_Tsu_text num2str(point.time_sec(k,1)) ',' num2str(point.T_ptc_su(k,1)) '; '];
        point.time_vs_Psu_text = [point.time_vs_Psu_text num2str(point.time_sec(k,1)) ',' num2str(point.P_ptc_su(k,1)) '; '];
        
    end
    
    eval([ Point_name ' = point;'])
    eval(['save(''' Point_name '.mat'', ''' Point_name ''' ) '])
end