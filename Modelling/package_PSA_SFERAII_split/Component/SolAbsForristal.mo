within package_PSA_SFERAII_split.Component;
model SolAbsForristal
  " PSA_SFERAII - 1D radial energy balance around the Heat Transfer Element based on the Forristal"
// It solves the 1D radial energy balance around the Heat Transfer Element based on the Forristal model see Forristal NREL 2003.

/*************************************************************/
/************************* 1. INPUTS *************************/
/*************************************************************/

  Modelica.Blocks.Interfaces.RealInput DNI
    annotation (Placement(transformation(extent={{-120,-40},{-80,0}}),
        iconTransformation(extent={{-106,-48},{-66,-8}})));
    Modelica.Blocks.Interfaces.RealInput v_wind
    annotation (Placement(transformation(extent={{-120,50},{-80,90}}),
        iconTransformation(extent={{-106,92},{-66,132}})));
  Modelica.Blocks.Interfaces.RealInput Theta "In radiants"
    annotation (Placement(transformation(extent={{-120,20},{-80,60}}),
        iconTransformation(extent={{-106,46},{-66,86}})));
  Modelica.Blocks.Interfaces.RealInput Tamb
    annotation (Placement(transformation(extent={{-120,-10},{-80,30}}),
        iconTransformation(extent={{-106,0},{-66,40}})));

/*************************************************************/
/*********************** 2. CONSTANTS  ***********************/
/*************************************************************/

constant Real pi = Modelica.Constants.pi "pi number";
constant Real Sigma = Modelica.Constants.sigma "Stefan-Boltzmann constant";
constant Real gg = Modelica.Constants.g_n
    "Standard acceleration of gravity on earth";

/*************************************************************/
/*********************** 3. PARAMETER ************************/
/*************************************************************/

/**** 3.1 Optical Properties ***********************/
parameter Real eps1 "HCE Shadowing [-]" annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real eps2 "Tracking error[-] "
                                        annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real eps3 "Geometry error [-]"
                                        annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real rho_cl "Mirror reflectivity [-]"
                                               annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real eps4 "Dirt on Mirrors [-] "
                                          annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real eps5 "Dirt on HCE [-]" annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real eps6 "Unaccounted [-]"
                                     annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real Tau_g = 0.963 "Glass Transmissivity [-]" annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real Alpha_g = 0.02 "Glass Absorptivity [-] " annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real Eps_g = 0.89 "Glass emissivity [-] " annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real Alpha_t =  0.97 "Tube Absorptivity [-]"
                                                      annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real a1_IAM = 0
    "IAM coefficiencs : IAM = 1-(a1_IAM*theta+a2_IAM*theta^2)/cos_theta"                       annotation (Dialog(group="Optical Properties", tab="General"));
parameter Real a2_IAM = 0
    "IAM coefficiencs : IAM = 1-(a1_IAM*theta+a2_IAM*theta^2)/cos_theta"                       annotation (Dialog(group="Optical Properties", tab="General"));
Real Eps_t[N] "Coating emissivity [-]";

/**** 3.2 PTC Properties ***********************/
parameter Integer N = 2 "number of cells";
parameter Modelica.SIunits.Length L "length of one tubes [m]";
parameter Modelica.SIunits.Length A_P "Aperture of the parabola [m]";

/**** 3.3 Envelope properties***/
parameter Modelica.SIunits.Length Dext_g = 0.12 "External glass diameter [m]" annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
parameter Modelica.SIunits.Length th_g = 0.0025 "Glass wall thickness [m]" annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
final parameter Modelica.SIunits.Length rext_g = Dext_g/2
    "Out radius glass [m]"                                                       annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
final parameter Modelica.SIunits.Length rint_g = rext_g-th_g
    "Int rad glass [m]"                                                          annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
final parameter Modelica.SIunits.Length Dint_g = 2*rint_g "Int rad glass [m]" annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
parameter Boolean GlassUD = false
    "if true, constant properties defined by the user"                                annotation (Dialog(group="Properties of the glass envelope", tab="General"));
parameter Modelica.SIunits.Density rho_g_ud = 2400 "Glass density [kg/m3]"  annotation (Dialog(enable=GlassUD,group="Properties of the glass envelope", tab="General"));
parameter Modelica.SIunits.SpecificHeatCapacity Cp_g_ud = 753
    "Specific heat capacity of the glass [J/kg.K]"                                                            annotation (Dialog(enable=GlassUD,group="Properties of the glass envelope", tab="General"));
parameter Modelica.SIunits.ThermalConductivity lambda_g_ud = 1.05
    "Thermal conductivity of the glass [W/m.K]"                                                               annotation (Dialog(enable=GlassUD,group="Properties of the glass envelope", tab="General"));
replaceable parameter Material.GlassEnvelope.MaterialBase  GlassMaterial constrainedby
    package_PSA_SFERAII_split.Component.Material.GlassEnvelope.MaterialBase
    "Select glass material for the envelope if not user defined"                                    annotation (choicesAllMatching=true, Dialog(enable=not
                                                                                                (GlassUD),group="Properties of the glass envelope", tab="General"));

/**** 3.4 Tube properties***/
parameter Modelica.SIunits.Length Dext_t = 0.07 "External diameter tube [m]" annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
parameter Modelica.SIunits.Length th_t= 0.002 "tube thickness [m]" annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
final parameter Modelica.SIunits.Length rext_t = Dext_t/2
    " Tube External Radius [m]"                                                          annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
final parameter Modelica.SIunits.Length rint_t= rext_t-th_t
    "Tube Internal Radius [m]"                                                           annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
parameter Boolean TubeUD = false
    "if true, constant properties defined by the user"                               annotation (Dialog(group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.Density rho_t_ud = 8000 "tube density [kg/m3]"  annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.SpecificHeatCapacity Cp_t_ud = 500
    "Specific heat capacity of the tube [J/kg.K] "                                                            annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.ThermalConductivity lambda_t_ud = 54
    "Thermal conductivity of the tube [W/m.K] "                                                             annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
replaceable parameter Material.TubeReceiver.MaterialBase TubeMaterial constrainedby
    package_PSA_SFERAII_split.Component.Material.TubeReceiver.MaterialBase
    "Select tube material for the receiver if not user defined"                                     annotation (choicesAllMatching=true, Dialog(enable=not
                                                                                                (TubeUD),group="Properties of the tube receiver", tab="General"));

/**** 3.5 Atmospheric properties***/
parameter Modelica.SIunits.Pressure Patm "Atmospheric pressure" annotation (Dialog(group="Atmospheric characteristic", tab="General"));
parameter Modelica.SIunits.ThermalConductivity k_air = 0.025574811
    "Thermal conductivity of air at constant pressure [W/m.K]"                                                                annotation (Dialog(group="Atmospheric characteristic", tab="General"));
parameter Modelica.SIunits.Density rho_air = 1.183650626
    "Density of air [kg/m3]"                                                      annotation (Dialog(group="Atmospheric characteristic", tab="General"));
parameter Modelica.SIunits.DynamicViscosity mu_air = 1.83055E-05
    "Dynamic viscosity of air at glass temperature"                                                               annotation (Dialog(group="Atmospheric characteristic", tab="General"));

Real C = if Re < 40 then
  0.75
 else
 if  Re > 40 and Re <1000 then
  0.51
 else
 if Re > 1000 and  Re < 200000 then
 0.26
 else
  0.076 "Coefficient for Nusselt number from Forristal ";
Real m = if Re < 40 then
  0.4
  else
  if  Re > 40 and Re <1000 then
  0.5
  else
  if Re > 1000 and  Re < 200000 then
  0.6
  else
  0.7 "Coefficient for Nusselt number from Forristal ";

Real n = if Pr < 10 then
     0.37
     else
     0.36 "Coefficient for Nusselt number from Forristal";
parameter Real Pr annotation (Dialog(group="Atmospheric characteristic", tab="General"));
Real Re "Reynolds number of atmosphere";
Real Nu "Nusselt number of atmophere";

/**** 3.6 Vaccum properties***/
parameter Modelica.SIunits.Pressure Pvacuum "Vacuum Pressure [Pa]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
final parameter Real P_mmhg = Pvacuum/133.322368 "Vacuum Pressure [mmHg]";
parameter Real GAMMA= 1.39 "Ratio of specific heats for the annulus gas [-]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
parameter Real DELTA = 3.53e-8 "Molecular diameter for the annulus gas [cm]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
parameter Real BB = 1.571 "Interaction coefficient [-]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
parameter Modelica.SIunits.ThermalConductivity k_st = 0.02551
    "Thermal conductivity at standard pressure and temperature"                                                            annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));

/*************************************************************/
/*********************** 4. Initial conditions****************/
/*************************************************************/

/**** 4.1 Glass temperature***/
parameter Modelica.SIunits.Temperature T_g_start_in
    "Temperature of the glass_inlet"                                                   annotation(Dialog(tab = "Initialization"));
parameter Modelica.SIunits.Temperature T_g_start_out
    "Temperature of the glass outlet"                                                  annotation(Dialog(tab = "Initialization"));
parameter Modelica.SIunits.Temperature T_g_start[N] = linspace(T_g_start_in,T_g_start_out,N)
    "start value of the tube temperature vector"                                                                                           annotation(Dialog(tab = "Initialization"));

/**** 4.2 Tube temperature***/
parameter Modelica.SIunits.Temperature T_t_start_in
    "Temperature of the tube inlet"                                                   annotation(Dialog(tab = "Initialization"));
parameter Modelica.SIunits.Temperature T_t_start_out
    "Average temperature of the tube outlet"                                                   annotation(Dialog(tab = "Initialization"));
parameter Modelica.SIunits.Temperature T_t_start[N] = linspace(T_t_start_in,T_t_start_out,N)
    "start value of the tube temperature vector"                                                                                           annotation(Dialog(tab = "Initialization"));

/*************************************************************/
/*********************** 5. VARIABLES ************************/
/*************************************************************/

/**** 5.1 Area of Collector and Reflector ***/
Modelica.SIunits.Area A_ext_g "Lateral external area of the glass";
Modelica.SIunits.Area A_ext_t "Lateral external area of the tube";
Modelica.SIunits.Area A_ref "Area of the reflector";
Modelica.SIunits.Area Am_t "Area of the metal cross section";
Modelica.SIunits.Area Am_g "Area of the glass cross section";

/**** 5.2 Heat transfer coefficients ***/
Modelica.SIunits.CoefficientOfHeatTransfer Gamma_vacuum[N]
    "Coefficient of heat transfer - vacuum";
Modelica.SIunits.CoefficientOfHeatTransfer Gamma_air
    "coefficient of heat transfer - ambient air";
Real LAMBDA[N] "mean free-path between collisions of a molecule [m]";

/**** 5.3 Temperature ***/
Modelica.SIunits.Temperature Tsky "Sky temperature";
Modelica.SIunits.Temperature T_int_t[N] "temperature of the internal tube";
Modelica.SIunits.Temperature T_t[N](start = T_t_start)
    "temperature at the center of the tube";
Modelica.SIunits.Temperature T_ext_t[N] "temperature of the external tube";
Modelica.SIunits.Temperature T_g_t[N]
    "Average temperature of the inner glass and the external tube";
Modelica.SIunits.Temperature T_int_g[N] "temperature of the internal glass";
Modelica.SIunits.Temperature T_g[N](start = T_g_start)
    "temperature at the center of the glass";
Modelica.SIunits.Temperature T_ext_g[N] "temperature of the external glass";

/**** 5.4 Heat power ***/
Modelica.SIunits.HeatFlowRate Q_sol_tot;
Modelica.SIunits.HeatFlowRate Q_glass_tot;
Modelica.SIunits.HeatFlowRate Q_tube_tot;
Modelica.SIunits.HeatFlowRate Q_loss_tot;
Modelica.SIunits.HeatFlowRate Q_abs[N];
Modelica.SIunits.HeatFlowRate Q_abs_tot;

/**** 5.5 Heat flux ***/
Modelica.SIunits.HeatFlux Phi_glass_tot "heat flux absorbed by the glass";
Modelica.SIunits.HeatFlux Phi_glass_tot_N[N]
    "Heat flux from the sun to the glass at each node";
Modelica.SIunits.HeatFlux Phi_tube_tot "Heat flux absorbed by the tube";
Modelica.SIunits.HeatFlux Phi_tube_tot_N[N]
    "Heat flux from the sun to the tube at each node";
Modelica.SIunits.HeatFlux Phi_tube_int[N] "Heat flux of tube internal";
Modelica.SIunits.HeatFlux Phi_tube_ext[N] "Heat flux of the tube external";
Modelica.SIunits.HeatFlux Phi_conv_gas[N] "heat flux convection in the vacuum";
Modelica.SIunits.HeatFlux Phi_rad_gas[N] "heat flux radiation in the vacuum";
Modelica.SIunits.HeatFlux Phi_glass_int[N] "heat flux of the glass internal";
Modelica.SIunits.HeatFlux Phi_glass_ext[N] "heat flux of the glass external";
Modelica.SIunits.HeatFlux Phi_conv_air[N] "Heat flux convection in the air";
Modelica.SIunits.HeatFlux Phi_rad_air[N] "heat flux radiation in the air";
Modelica.SIunits.HeatFlux Phi_loss
    "Heat losses with respect to the reflector surface";

/**** 5.6 Efficiencies and other variables***/
Real eta_opt;
Real eta_opt_g;
Real eta_opt_t;
Real IAM "Incident Angle Modifier";
Real eta_th[N] "Thermal efficiency at each node";
Real eta_tot[N] "Total efficiency : eta_th * eta_opt_t";
Real Eta_th "Average Thermal efficiency";
Real Eta_tot "Average Total efficiency";
Real cos_theta "cos theta";

/**** 5.7 Materials properties ***/
Modelica.SIunits.Density rho_g[N];
Modelica.SIunits.SpecificHeatCapacity Cp_g[N];
Modelica.SIunits.ThermalConductivity lambda_g[N];
Modelica.SIunits.Density rho_t[N];
Modelica.SIunits.SpecificHeatCapacity Cp_t[N];
Modelica.SIunits.ThermalConductivity lambda_t[N];

/**** 5.8 Thermal port ***/
 ThermoCycle.Interfaces.HeatTransfer.ThermalPort wall_int(N=N) annotation (
      Placement(transformation(extent={{60,20},{80,40}}), iconTransformation(
          extent={{140,14},{160,34}})));

/*************************************************************/
/************************* 6. MODELLING***********************/
/*************************************************************/
  Modelica.Blocks.Interfaces.IntegerInput Focus(start=1)
    "PTC focused if Focus = 1, PTC unfocused = 0" annotation (Placement(
        transformation(extent={{-106,-96},{-66,-56}}), iconTransformation(
          extent={{-106,-96},{-66,-56}})));
equation

//Sky temperature //
Tsky = Tamb - 8 "Assumption from Forristal";

//Cos_theta and Incidence angle modifier //
cos_theta = Modelica.Math.cos(Theta);
IAM = 1-(a1_IAM*Theta+a2_IAM*Theta^2)/cos_theta
    "Correlation based on Vazuela's paper";

//Optical efficiency //
eta_opt = eps1*eps2*eps3*eps4*eps5*eps6*rho_cl*IAM;

//Optical efficiency glass //
eta_opt_g = eta_opt * Alpha_g;

//Optical efficiency tube //
eta_opt_t = eta_opt * Alpha_t * Tau_g;

//Total area of the reflector //
A_ref = (Focus*L*A_P) + (1-Focus)*L*Dext_t;

//Total area of the external glass//
A_ext_g = pi*Dext_g*L;

//Total area of the external tube //
A_ext_t = pi*Dext_t*L;

//Total thermal enery incident to te collector //
Q_sol_tot = DNI*cos_theta*A_ref;

//Total thermal energy on the glass from the Sun //
Q_glass_tot = eta_opt_g*DNI*cos_theta*A_ref
    "Thermal power concentrated on the glass envelope [W]";
Phi_glass_tot = Q_glass_tot / A_ext_g
    "Heat flux on the glass envelope [W/m2_g_ext]";

// Total thermal energy on the tube from the Sun//
Q_tube_tot = eta_opt_t*DNI*cos_theta*A_ref
    "Thermal power concentrated on the receiver tube [W]";
Phi_tube_tot = Q_tube_tot/A_ext_t "Heat flux on the receiver tube [W/m2_t_ext]";

//Conduction glass&tube datas//
assert(rext_t > rint_t, "External Radius must be greater than internal radius");
Am_t = (rext_t^2 - rint_t^2)*pi "Area of the metal cross section";
assert(rext_g > rint_g, "External Radius must be greater than internal radius");
Am_g = (rext_g^2 - rint_g^2)*pi "Area of the glass cross section";

//Convective Heat Transfer Coefficient air //
Re = v_wind*rho_air*Dext_g/mu_air "Reynolds of ambiant air";
Nu = C* Re^m * Pr ^ n "Nusselt of ambient air";
Gamma_air = k_air*Nu/Dext_g
    "Convective heat transfer coefficient of ambient air";

for i in 1:N loop

  //Tube Emissivity //
  Eps_t[i] = 0.06282 + (1.208E-4)*(T_ext_t[i]-273.15)+ (1.907E-7)*(T_ext_t[i]-273.15)^2
      "Corrected UVAC emissivity";

  //Heat transfer in vaccum  //
  T_g_t[i] = (T_int_g[i] + T_ext_t[i])/2 "mean temperature inside annulus";
  LAMBDA[i] = 2.331E-20*(T_g_t[i])/(P_mmhg * DELTA^2)
      "ratio of specific heat in annulus";
  Gamma_vacuum[i] = k_st /((Dext_t/2)*log(rint_g/rext_t)+BB*LAMBDA[i]*(rext_t/rint_g+1));

  //Heat flux on tube //
  Phi_tube_tot_N[i] = Phi_tube_tot;

  //Heat flux on glass //
  Phi_glass_tot_N[i] = Phi_glass_tot;

  //Convection to the external air//
  Phi_conv_air[i] = Gamma_air * (T_ext_g[i] - Tamb);

  //Radiation to the external air//
  Phi_rad_air[i] = Eps_g * Sigma * (T_ext_g[i]^4 - Tsky^4);

  //Heat balance on the glass external envelop//
  Phi_glass_ext[i] = Phi_glass_tot_N[i] - Phi_rad_air[i] - Phi_conv_air[i];

  //Conduction in the glass //
  if GlassUD then
    rho_g[i] = rho_g_ud;
    Cp_g[i] = Cp_g_ud;
    lambda_g[i] = lambda_g_ud;
 else
    rho_g[i] = GlassMaterial.a0_rho_g + GlassMaterial.a1_rho_g*(T_g[i]-273.15) + GlassMaterial.a2_rho_g*(T_g[i]-273.15) ^2;
    Cp_g[i] = GlassMaterial.a0_cp_g + GlassMaterial.a1_cp_g*(T_g[i]-273.15) + GlassMaterial.a2_cp_g*(T_g[i]-273.15) ^2;
    lambda_g[i] = GlassMaterial.a0_cp_g + GlassMaterial.a1_lambda_g*(T_g[i]-273.15) + GlassMaterial.a2_lambda_g*(T_g[i]-273.15) ^2;
  end if;
  rho_g[i]*Cp_g[i]*Am_g*der(T_g[i]) = rint_g*2*pi*Phi_glass_int[i] + rext_g*2*pi*Phi_glass_ext[i]
      "Linear energy balance i.e. [W/m_glass]";
  Phi_glass_ext[i] = lambda_g[i]/(rext_g*log((2*rext_g)/(rint_g + rext_g)))*(T_ext_g[i] - T_g[i])
      "Heat conduction through the external half-thickness";
  Phi_glass_int[i] = lambda_g[i]/(rint_g*log((rint_g + rext_g)/(2*rint_g)))*(T_int_g[i] - T_g[i])
      "Heat conduction through the internal half-thickness";

  //Connection Internal Heat flow to the glass //
  Phi_glass_int[i] = (Phi_conv_gas[i] + Phi_rad_gas[i]);

  //Convection in the vacuum //
  Phi_conv_gas[i] = Gamma_vacuum[i] *(T_ext_t[i] - T_int_g[i]);

  //Radiation in the vacuum //
  Phi_rad_gas[i] = Sigma*(T_ext_t[i]^4 - T_int_g[i]^4)/(1/Eps_t[N] + Dext_t/Dint_g*(1/Eps_g-1));

  //Heat flux to the tube //
  Phi_tube_ext[i] = Phi_tube_tot_N[i] - Phi_glass_int[i];

  //Conduction tube //
  if TubeUD then
    rho_t[i] = rho_t_ud;
    Cp_t[i] = Cp_t_ud;
    lambda_t[i] = lambda_t_ud;
  else
    rho_t[i] = TubeMaterial.a0_rho_t + TubeMaterial.a1_rho_t*(T_t[i]-273.15) + TubeMaterial.a2_rho_t*(T_t[i]-273.15) ^2;
    Cp_t[i] = TubeMaterial.a0_cp_t + TubeMaterial.a1_cp_t*(T_t[i]-273.15) + TubeMaterial.a2_cp_t*(T_t[i]-273.15) ^2;
    lambda_t[i] = TubeMaterial.a0_cp_t + TubeMaterial.a1_lambda_t*(T_t[i]-273.15) + TubeMaterial.a2_lambda_t*(T_t[i]-273.15) ^2;
  end if;
  rho_t[i]*Cp_t[i]*Am_t*der(T_t[i]) = rint_t*2*pi*Phi_tube_int[i] + rext_t*2*pi*Phi_tube_ext[i]
      "Energy balance [W/m_tube]";
  Phi_tube_ext[i] = lambda_t[i]/(rext_t*log((2*rext_t)/(rint_t + rext_t)))*(T_ext_t[i] - T_t[i])
      "Heat conduction through the external half-thickness";
  Phi_tube_int[i] = lambda_t[i]/(rint_t*log((rint_t + rext_t)/(2*rint_t)))*(T_int_t[i] - T_t[i])
      "Heat conduction through the internal half-thickness";

  //Fluid interaction //
  wall_int.phi[i]=Phi_tube_int[i];
  T_int_t[i] = wall_int.T[i];

  // THERMAL EFFICIENCY AT EACH NODE //
  Q_abs[i] = Phi_tube_int[i]*2*rint_t*pi*L/N
      "Heat power absorbed by the fluid in each cell";
  if Q_tube_tot > 0 then
    eta_th[i] = Q_abs[i]/(Q_tube_tot/N);
    eta_tot[i] = eta_th[i] *eta_opt_t;
  else
    eta_th[i] = 0;
    eta_tot[i] = eta_th[i] * eta_opt_t;
  end if;

end for;

// THERMAL LOSSES PER REFLECTOR SURFACE [W/m2]
Q_loss_tot = (sum(Phi_rad_air) + sum(Phi_conv_air))*(A_ext_g /N) "[W]";
Q_abs_tot = sum(Q_abs[:]) "[W]";
Phi_loss =  (sum(Phi_rad_air) + sum(Phi_conv_air))*A_ext_g /(A_ref*N) "[W/m2]";

//THERMAL EFFICIENCY
Eta_th = sum(eta_th)/N;

//TOTAL EFFICIENCY
Eta_tot = Eta_th*eta_opt_t;
                                                                                                      annotation(Dialog(tab = "Initialisation"),
             Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{160,140}}),
                     graphics),
    Icon(coordinateSystem(extent={{-100,-100},{160,140}}, preserveAspectRatio=false),
         graphics),Documentation(info="<HTML>
          
         <p><big>Model <b>SolAbsForristal</b>  represents the one-dimensional radial energy balance between the Heat Collector Element (HCE) and the atmosphere based on the
          <a href=\"http://www.nrel.gov/csp/troughnet/pdfs/34169.pdf\">Forristal model</a>.
      <p><big> The terms in the energy balance depends on the collector type, the HCE condition, the optical properties and the ambient condition.
         
         <p><big>The phenomena represented by the model are:
         <ul><li>Convection in the heat transfer fluid.
         <li>Conduction and thermal energy storage in the metal pipe.
         <li> Convection and radiation transfer in the vacuum between the glass envelope and the metal pipe.
         <li> Conduction and thermal energy storage in the glass envelope. 
         <li> Radiation and convection to the environment.
         </ul>
         
         <p><big>The model assumptions are:
<ul><li> Temperatures, thermal energy flux and thermodynamic properties are considered uniform around the circumference of the HCE.
 <li> Solar absorption is treated as a linear phenomenon.
 <li> Constant density in the metal pipe and in the glass envelope.
 <li> Constant heat capacity in the metal pipe and in the glass envelope.          
        </ul> 
         
</HTML>"));
end SolAbsForristal;
