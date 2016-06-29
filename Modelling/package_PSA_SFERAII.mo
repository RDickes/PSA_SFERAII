within ;
package package_PSA_SFERAII
  "package including all the models required for the PSA_SFERAII_project"
  package Simulations "folder including the different simulation"
    model Basic1 "First model - most basic"
      ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot(
          redeclare package Medium =
            ThermoCycle.Media.Incompressible.IncompressibleCP.HighTemperature.Therminol_66,
          Mdot_0=1)
        annotation (Placement(transformation(extent={{-82,-82},{-62,-62}})));
      ThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP(redeclare package
          Medium =
            ThermoCycle.Media.Incompressible.IncompressibleCP.HighTemperature.Therminol_66,
          p0=200000)
        annotation (Placement(transformation(extent={{60,70},{80,90}})));
      ThermoCycle.Components.Units.Solar.SolarField_Forristal_Inc
        solarField_Forristal_Inc(
        redeclare package Medium1 =
            ThermoCycle.Media.Incompressible.IncompressibleCP.HighTemperature.Therminol_66,
        T_g_start_in=373.15,
        T_g_start_out=373.15,
        T_t_start_in=373.15,
        T_t_start_out=373.15,
        Tstart_inlet=373.15,
        Tstart_outlet=373.15,
        pstart=200000,
        Discretization=ThermoCycle.Functions.Enumerations.Discretizations.upwind_AllowFlowReversal)
        annotation (Placement(transformation(extent={{-38,-28},{22,38}})));

      Modelica.Blocks.Sources.Step v_wind_input(height=0)
        annotation (Placement(transformation(extent={{-90,66},{-76,80}})));
      Modelica.Blocks.Sources.Step theta_input(height=0)
        annotation (Placement(transformation(extent={{-92,34},{-76,50}})));
      Modelica.Blocks.Sources.Step T_amb_input(height=298)
        annotation (Placement(transformation(extent={{-92,2},{-78,16}})));
      Modelica.Blocks.Sources.Step DNI_input(height=1000)
        annotation (Placement(transformation(extent={{-92,-28},{-78,-14}})));
    equation
      connect(v_wind_input.y, solarField_Forristal_Inc.v_wind) annotation (Line(
          points={{-75.3,73},{-55.65,73},{-55.65,36.02},{-32.6667,36.02}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(theta_input.y, solarField_Forristal_Inc.Theta) annotation (Line(
          points={{-75.2,42},{-62,42},{-62,15.56},{-32,15.56}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_amb_input.y, solarField_Forristal_Inc.Tamb) annotation (Line(
          points={{-77.3,9},{-57.65,9},{-57.65,-1.6},{-32.6667,-1.6}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(DNI_input.y, solarField_Forristal_Inc.DNI) annotation (Line(
          points={{-77.3,-21},{-57.65,-21},{-57.65,-18.76},{-32.6667,-18.76}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sourceMdot.flangeB, solarField_Forristal_Inc.InFlow) annotation (
          Line(
          points={{-63,-72},{2,-72},{2,-28},{1.33333,-28}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(solarField_Forristal_Inc.OutFlow, sinkP.flangeB) annotation (Line(
          points={{2,43.28},{2,80},{61.6,80}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics));
    end Basic1;
  end Simulations;

  package Media "Include new media data"
    package Sytherm800 "Sytherm800 Incompressible - TableBased"
      extends ThermoCycle.Media.Incompressible.TableBased(
        mediumName="Syltherm800",
        T_min=Modelica.SIunits.Conversions.from_degC(-40),
        T_max=Modelica.SIunits.Conversions.from_degC(400),
        TinK=false,
        T0=273.15,
        tableDensity = [-40, 990.61; 0, 953.16; 40, 917.07; 80, 881.68; 120, 846.35; 160, 810.45; 200, 773.33; 240, 734.35; 280, 692.87; 320, 648.24; 360, 599.83; 400, 547],
        tableHeatCapacity = [-40, 1506; 0, 1574; 40, 1643; 80, 1711; 120, 1779; 160, 1847; 200, 1916; 240, 1984; 280, 2052; 320, 2121; 360, 2189; 400, 2257],
        tableConductivity = [-40, 0.1463; 0, 0.1388; 40, 0.1312; 80, 0.1237; 120, 0.1162; 160, 0.1087; 200, 0.1012; 240, 0.0936; 280, 0.0861; 320, 0.0786; 360, 0.0711; 400, 0.0635],
        tableViscosity = [-40, 0.05105; 0, 0.01533; 40, 0.007; 80, 0.00386; 120, 0.00236; 160, 0.00154; 200, 0.00105; 240, 0.00074; 280, 0.00054; 320, 0.00041; 360, 0.00031; 400, 0.00025],
        tableVaporPressure = [-40, 0; 0, 0; 40, 100; 80, 1460; 120, 9300; 160, 35000; 200, 94600; 240, 204800; 280, 380200; 320, 630500; 360, 961200; 400, 1373000]);

       // Density ---->  [kg/m3]
       // HeatCapacity ----> [J/kg/K]
       // Conuctivity  ----> [W/m/K]
       // Viscosity  ---->    [Pa.s]
       // Vapor pressure ---->  [Pa]
    end Sytherm800;
  end Media;

  package Component
    model SolAbsForristal
      " PSA_SFERAII - 1D radial energy balance around the Heat Transfer Element based on the Forristal"
    // It solves the 1D radial energy balance around the Heat Transfer Element based on the Forristal model see Forristal NREL 2003.

    /*************************************************************/
    /************************* 1. INPUTS *************************/
    /*************************************************************/

      Modelica.Blocks.Interfaces.RealInput DNI
        annotation (Placement(transformation(extent={{-120,-40},{-80,0}}),
            iconTransformation(extent={{-108,-96},{-68,-56}})));
        Modelica.Blocks.Interfaces.RealInput v_wind
        annotation (Placement(transformation(extent={{-120,50},{-80,90}}),
            iconTransformation(extent={{-106,64},{-66,104}})));
      Modelica.Blocks.Interfaces.RealInput Theta "In radiants"
        annotation (Placement(transformation(extent={{-120,20},{-80,60}}),
            iconTransformation(extent={{-106,2},{-66,42}})));
      Modelica.Blocks.Interfaces.RealInput Tamb
        annotation (Placement(transformation(extent={{-120,-10},{-80,30}}),
            iconTransformation(extent={{-106,-46},{-66,-6}})));

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
    parameter Modelica.SIunits.Length Dext_g = 0.12
        "External glass diameter [m]"                                             annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.Length th_g = 0.0025 "Glass wall thickness [m]" annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    final parameter Modelica.SIunits.Length rext_g = Dext_g/2
        "Out radius glass [m]"                                                       annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    final parameter Modelica.SIunits.Length rint_g = rext_g-th_g
        "Int rad glass [m]"                                                          annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.Density rho_g "Glass density [kg/m3] " annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.SpecificHeatCapacity Cp_g
        "Specific heat capacity of the glass [J/kg.K]"                                                  annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity lambda_g
        "Thermal conductivity of the glass [W/m.K]"                                                     annotation (Dialog(group="GeometriesAndProperties of the glass envelope", tab="General"));

    /**** 3.4 Tube properties***/
    parameter Modelica.SIunits.Length Dext_t = 0.07
        "External diameter tube [m]"                                             annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    parameter Modelica.SIunits.Length th_t= 0.002 "tube thickness [m]" annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    final parameter Modelica.SIunits.Length rext_t = Dext_t/2
        " Tube External Radius [m]"                                                          annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    final parameter Modelica.SIunits.Length rint_t= rext_t-th_t
        "Tube Internal Radius [m]"                                                           annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    parameter Modelica.SIunits.Density rho_t "tube density [kg/m3]" annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    parameter Modelica.SIunits.SpecificHeatCapacity Cp_t
        "Specific heat capacity of the tube [J/kg.K]"                                                    annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity lambda_t
        "Thermal conductivity of the tube [W/m.K]"                                                      annotation (Dialog(group="GeometriesAndProperties of the metal envelope", tab="General"));

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
    parameter Real GAMMA= 1.39
        "Ratio of specific heats for the annulus gas [-]"                        annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
    parameter Real DELTA = 3.53e-8
        "Molecular diameter for the annulus gas [cm]"                            annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
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
    Modelica.SIunits.HeatFlowRate Q_glass_tot;
    Modelica.SIunits.HeatFlowRate Q_tube_tot;
    Modelica.SIunits.HeatFlowRate Q_abs[N];

    /**** 5.5 Heat flux ***/
    Modelica.SIunits.HeatFlux Phi_glass_tot "heat flux absorbed by the glass";
    Modelica.SIunits.HeatFlux Phi_glass_tot_N[N]
        "Heat flux from the sun to the glass at each node";
    Modelica.SIunits.HeatFlux Phi_tube_tot "Heat flux absorbed by the tube";
    Modelica.SIunits.HeatFlux Phi_tube_tot_N[N]
        "Heat flux from the sun to the tube at each node";
    Modelica.SIunits.HeatFlux Phi_tube_int[N] "Heat flux of tube internal";
    Modelica.SIunits.HeatFlux Phi_tube_ext[N] "Heat flux of the tube external";
    Modelica.SIunits.HeatFlux Phi_conv_gas[N]
        "heat flux convection in the vacuum";
    Modelica.SIunits.HeatFlux Phi_rad_gas[N]
        "heat flux radiation in the vacuum";
    Modelica.SIunits.HeatFlux Phi_glass_int[N]
        "heat flux of the glass internal";
    Modelica.SIunits.HeatFlux Phi_glass_ext[N]
        "heat flux of the glass external";
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
    Real eta_TOT[N] "Total efficiency : eta_th * eta_opt_t";
    Real Eta_th "Average Thermal efficiency";
    Real Eta_TOT "Average Total efficiency";
    Real cos_theta "cos theta";

    /**** 5.7 Thermal port ***/
     ThermoCycle.Interfaces.HeatTransfer.ThermalPort wall_int(N=N) annotation (
          Placement(transformation(extent={{60,20},{80,40}}), iconTransformation(
              extent={{80,-10},{100,10}})));

    /*************************************************************/
    /************************* 6. MODELLING***********************/
    /*************************************************************/
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
    A_ref = L*A_P;

    //Total area of the external glass//
    A_ext_g = pi*Dext_g*L;

    //Total area of the external tube //
    A_ext_t = pi*Dext_t*L;

    //Total thermal energy on the glass from the Sun //
    Q_glass_tot = eta_opt_g*DNI*cos_theta*A_ref
        "Thermal power concentrated on the glass envelope [W]";
    Phi_glass_tot = Q_glass_tot / A_ext_g
        "Heat flux on the glass envelope [W/m2_g_ext]";

    // Total thermal energy on the tube from the Sun//
    Q_tube_tot = eta_opt_t*DNI*cos_theta*A_ref
        "Thermal power concentrated on the receiver tube [W]";
    Phi_tube_tot = Q_tube_tot/A_ext_t
        "Heat flux on the receiver tube [W/m2_t_ext]";

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
      Eps_t[i] = 0.062 + (2E-7)*(T_ext_t[i]-273.15)^2 "--> TO CHANGE";

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
      rho_g*Cp_g*Am_g*der(T_g[i]) = rint_g*2*pi*Phi_glass_int[i] + rext_g*2*pi*Phi_glass_ext[i]
          "Linear energy balance i.e. [W/m_glass]";
      Phi_glass_ext[i] = lambda_g/(rext_g*log((2*rext_g)/(rint_g + rext_g)))*(T_ext_g[i] - T_g[i])
          "Heat conduction through the external half-thickness";
      Phi_glass_int[i] = lambda_g/(rint_g*log((rint_g + rext_g)/(2*rint_g)))*(T_int_g[i] - T_g[i])
          "Heat conduction through the internal half-thickness";

      //Connection Internal Heat flow to the glass //
      Phi_glass_int[i] = (Phi_conv_gas[i] + Phi_rad_gas[i]);

      //Convection in the vacuum //
      Phi_conv_gas[i] = Gamma_vacuum[i] *(T_ext_t[i] - T_int_g[i]);

      //Radiation in the vacuum //
      Phi_rad_gas[i] = Sigma*(T_ext_t[i]^4 - T_int_g[i]^4)/(1/Eps_t[N] + Dext_t/Dext_g*(1/Eps_g-1))
          "FAUTE --> c'est Dext_t/Dint_g";

      //Heat flux to the tube //
      Phi_tube_ext[i] = Phi_tube_tot_N[i] - Phi_glass_int[i];

      //Conduction tube //
      rho_t*Cp_t*Am_t*der(T_t[i]) = rint_t*2*pi*Phi_tube_int[i] + rext_t*2*pi*Phi_tube_ext[i]
          "Energy balance [W/m_tube]";
      Phi_tube_ext[i] = lambda_t/(rext_t*log((2*rext_t)/(rint_t + rext_t)))*(T_ext_t[i] - T_t[i])
          "Heat conduction through the external half-thickness";
      Phi_tube_int[i] = lambda_t/(rint_t*log((rint_t + rext_t)/(2*rint_t)))*(T_int_t[i] - T_t[i])
          "Heat conduction through the internal half-thickness";

      //Fluid interaction //
      wall_int.phi[i]=Phi_tube_int[i] "should not be a - ? --> to check";
      T_int_t[i] = wall_int.T[i];

      // THERMAL EFFICIENCY AT EACH NODE //
      Q_abs[i] = Phi_tube_int[i]*2*rint_t*pi*L/N
          "Heat power absorbed by the fluid in each cell";
      if Q_tube_tot > 0 then
        eta_th[i] = Q_abs[i]/(Q_tube_tot/N);
        eta_TOT[i] = eta_th[i] *eta_opt_t;
      else
        eta_th[i] = 0;
        eta_TOT[i] = eta_th[i] * eta_opt_t;
      end if;

    end for;

    // THERMAL LOSSES PER REFLECTOR SURFACE [W/m2]
    Phi_loss =  (sum(Phi_rad_air) + sum(Phi_conv_air))*A_ext_g /(A_ref*N);

    //THERMAL EFFICIENCY
    Eta_th = sum(eta_th)/N;

    //TOTAL EFFICIENCY
    Eta_TOT = Eta_th*eta_opt_t;
                                                                                                          annotation(Dialog(tab = "Initialisation"),
                 Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                         graphics),
        Icon(graphics),Documentation(info="<HTML>
          
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

    model SolarField_Forristal_Inc
      "PSA_SFERAII - Solar field model with collectors based on Forristal model for incompressible fluids"

    /*************************************************************/
    /************************* 1. FLUID **************************/
    /*************************************************************/
      replaceable package Medium1 = ThermoCycle.Media.DummyFluid constrainedby
        Modelica.Media.Interfaces.PartialMedium                                                                        annotation (choicesAllMatching = true);

    /*************************************************************/
    /*********************** 2. CONSTANTS  ***********************/
    /*************************************************************/

    constant Real  pi = Modelica.Constants.pi;

    /*************************************************************/
    /*********************** 3. PARAMETER ************************/
    /*************************************************************/

    /**** 3.1 Optical Properties ***********************/
    parameter Real eps1 = 0.9754 "HCE Shadowing [-]"
                                                    annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real eps2 = 0.994 "Tracking error [-]"
                                                    annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real eps3 = 0.98 "Geometry error [-]"
                                                   annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real eps4 = 0.962566845 "Dirt on Mirrors [-]"
                                                           annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real eps5 = 0.981283422 "Dirt on HCE [-]"
                                                       annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real eps6 = 0.96 "Unaccounted [-]"
                                                annotation (Dialog(group="Optical Properties", tab="General"));
    parameter Real rho_cl = 0.935 "Mirror reflectivity [-]"
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

    /**** 3.2 PTC Properties ***********************/
    parameter Integer N(min=1) = 10 "Number of cells per tube";
    parameter Integer Ns(min=1) = 2 "Number of tube in series";
    parameter Integer Nt(min=1) = 1 "Number of tubes in parallel";
    parameter Modelica.SIunits.Length L = 8 "Length of one tube [m]";
    parameter Modelica.SIunits.Length A_P = 5 "Aperture of the parabola [m]";

    final parameter Modelica.SIunits.Area A_lateral= L*D_int_t*pi
        "Lateral internal surface of the metal tube for one collector";
    final parameter Modelica.SIunits.Volume V_tube_int = pi*D_int_t^2/4*L
        "Internal volume of the metal tube for one collector";

    /**** 3.3 Envelope properties***/
    parameter Modelica.SIunits.Length Dext_g = 0.12
        "External glass diameter [m] "                                            annotation (Dialog(group="Properties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.Length th_g = 0.0025 "Glass thickness [m] "
                                                                          annotation (Dialog(group="Properties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.Density rho_g = 2400 "Glass density [kg/m3]"    annotation (Dialog(group="Properties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.SpecificHeatCapacity Cp_g = 753
        "Specific heat capacity of the glass [J/kg.K]"                                                         annotation (Dialog(group="Properties of the glass envelope", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity lambda_g = 1.05
        "Thermal conductivity of the glass [W/m.K]"                                                            annotation (Dialog(group="Properties of the glass envelope", tab="General"));

    /**** 3.4 Tube properties***/
    parameter Modelica.SIunits.Length Dext_t =  0.07
        "External diameter tube [m] "                                             annotation (Dialog(group="Properties of the tube receiver", tab="General"));
    parameter Modelica.SIunits.Length th_t =  0.002 "tube thickness [m] "
                                                                         annotation (Dialog(group="Properties of the tube receiver", tab="General"));
    final parameter Modelica.SIunits.Length D_int_t= Dext_t - 2*th_t
        "internal diameter of the metal tube [m] ";
    parameter Modelica.SIunits.Density rho_t = 8000 "tube density [kg/m3]"  annotation (Dialog(group="Properties of the tube receiver", tab="General"));
    parameter Modelica.SIunits.SpecificHeatCapacity Cp_t = 500
        "Specific heat capacity of the tube [J/kg.K] "                                                         annotation (Dialog(group="Properties of the tube receiver", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity lambda_t = 54
        "Thermal conductivity of the tube [W/m.K] "                                                          annotation (Dialog(group="Properties of the tube receiver", tab="General"));

    /**** 3.5 Atmospheric properties***/
    parameter Real Pr= 0.72 "Prandt number [-]"
                                               annotation (Dialog(group="Atmospheric characteristic", tab="General"));
    parameter Modelica.SIunits.Pressure Patm = 1e5 "Atmospheric Pressure [Pa]"   annotation (Dialog(group="Atmospheric characteristic", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity k_air = 0.025574811
        "Thermal conductivity of ambient air [W/m.K]"                                                                annotation (Dialog(group="Atmospheric characteristic", tab="General"));
    parameter Modelica.SIunits.Density rho_air = 1.183650626
        "Density of ambient air [kg/m3]"                                                      annotation (Dialog(group="Atmospheric characteristic", tab="General"));
    parameter Modelica.SIunits.DynamicViscosity mu_air = 1.83055E-05
        "Dynamic viscosity of ambient air"                                                               annotation (Dialog(group="Atmospheric characteristic", tab="General"));

    /**** 3.6 Vaccum properties***/
    parameter Modelica.SIunits.Pressure Pvacuum = 0.013332237
        "Vacuum Pressure [Pa]"                                                      annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
    parameter Real GAMMA= 1.39
        "Ratio of specific heats for the annulus gas [-]"                        annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
    parameter Real DELTA = 3.53e-8
        "Molecular diameter for the annulus gas [cm]"                            annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
    parameter Real BB = 1.571 "Interaction coefficient [-]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
    parameter Modelica.SIunits.ThermalConductivity k_st = 0.02551
        "Thermal conductivity at standard pressure and temperature"                                                            annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));

    /*************************************************************/
    /*********************** 4. Initial conditions****************/
    /*************************************************************/

    /**** 4.1 Glass temperature***/
    parameter Modelica.SIunits.Temperature T_g_start_in = 42 + 273.15
        "Temperature at the inlet of the glass"                                                                 annotation (Dialog(tab="Initialization"));
    parameter Modelica.SIunits.Temperature T_g_start_out = 54 + 273.15
        "Temperature at the outlet of the glass"                                                                 annotation (Dialog(tab="Initialization"));
    final parameter Modelica.SIunits.Temperature[Ns] T_g_start_inlet_collector =  ThermoCycle.Functions.Solar.T_start_inlet(T_start_inlet=T_g_start_in,T_start_outlet=T_g_start_out,Ns=Ns);
    final parameter Modelica.SIunits.Temperature[Ns] T_g_start_outlet_collector = ThermoCycle.Functions.Solar.T_start_outlet(T_start_inlet=T_g_start_in,T_start_outlet=T_g_start_out,Ns=Ns);

    /**** 4.2 Tube temperature***/
    parameter Modelica.SIunits.Temperature T_t_start_in = 195 + 273.15
        "Temperature at the inlet of the tube"                                                                  annotation (Dialog(tab="Initialization"));
    parameter Modelica.SIunits.Temperature T_t_start_out = 305 + 273.15
        "Temperature at the outlet of the tube"                                                                 annotation (Dialog(tab="Initialization"));
    final parameter Modelica.SIunits.Temperature[Ns] T_t_start_inlet_collector =  ThermoCycle.Functions.Solar.T_start_inlet(T_start_inlet=T_t_start_in,T_start_outlet=T_t_start_out,Ns=Ns);
    final parameter Modelica.SIunits.Temperature[Ns] T_t_start_outlet_collector = ThermoCycle.Functions.Solar.T_start_outlet(T_start_inlet=T_t_start_in,T_start_outlet=T_t_start_out,Ns=Ns);

    /**** 4.3 Fluid temperature***/
    parameter Modelica.SIunits.CoefficientOfHeatTransfer Unom=300
        "Coefficient of heat transfer"                                                            annotation (Dialog(group="Heat transfer", tab="General"));
    parameter Modelica.SIunits.MassFlowRate Mdotnom= 0.5
        "Total nominal Mass flow";
    parameter Modelica.SIunits.Temperature Tstart_inlet
        "Temperature of the fluid at the inlet of the collector"                                                  annotation (Dialog(tab="Initialization"));
    parameter Modelica.SIunits.Temperature Tstart_outlet
        "Temperature of the fluid at the outlet of the collector"                                                   annotation (Dialog(tab="Initialization"));
    final parameter Modelica.SIunits.Temperature[Ns] Tstart_inlet_collector =  ThermoCycle.Functions.Solar.T_start_inlet(T_start_inlet=Tstart_inlet,T_start_outlet=Tstart_outlet,Ns=Ns);
    final parameter Modelica.SIunits.Temperature[Ns] Tstart_outlet_collector = ThermoCycle.Functions.Solar.T_start_outlet(T_start_inlet=Tstart_inlet,T_start_outlet=Tstart_outlet,Ns=Ns);
    parameter Modelica.SIunits.Pressure pstart
        "Temperature of the fluid at the inlet of the collector"                                        annotation (Dialog(tab="Initialization"));

    /**** 4.4 Initial state : steady state  or not ***/
     parameter Boolean steadystate_T_fl=false
        "if true, sets the derivative of the fluid Temperature in each cell to zero during Initialization"
                                                                                                        annotation (Dialog(group="Initialization options", tab="Initialization"));

    /*************************************************************/
    /*********************** 5. Numerical options ****************/
    /*************************************************************/

    /**** 5.1 Discretization scheme ***/
    import ThermoCycle.Functions.Enumerations.Discretizations;
    parameter Discretizations Discretization=ThermoCycle.Functions.Enumerations.Discretizations.centr_diff
        "Selection of the spatial discretization scheme"                                                                                                     annotation (Dialog(tab="Numerical options"));

    /**** 5.2 Heat transfer model scheme ***/
    replaceable model FluidHeatTransferModel =
          ThermoCycle.Components.HeatFlow.HeatTransfer.MassFlowDependence
    constrainedby
        ThermoCycle.Components.HeatFlow.HeatTransfer.BaseClasses.PartialHeatTransferZones
        "Fluid heat transfer model"                                                                                             annotation (Dialog(group="Heat transfer", tab="General"),choicesAllMatching=true);

     /******************************* COMPONENTS ***********************************/

        package_PSA_SFERAII.Component.SolAbsForristal[Ns] solAbs(
        each eps1=eps1,
        each eps2=eps2,
        each eps3=eps3,
        each rho_cl=rho_cl,
        each eps4=eps4,
        each eps5=eps5,
        each eps6=eps6,
        each Tau_g = Tau_g,
        each Alpha_g = Alpha_g,
        each Eps_g = Eps_g,
        each Alpha_t = Alpha_t,
        each a1_IAM = a1_IAM,
        each a2_IAM = a2_IAM,
        each N=N,
        each L=L,
        each A_P=A_P,
        each rho_g=rho_g,
        each Cp_g=Cp_g,
        each lambda_g=lambda_g,
        each rho_t=rho_t,
        each Cp_t=Cp_t,
        each lambda_t=lambda_t,
        each Patm=Patm,
        each Pr=Pr,
        each k_air = k_air,
        each rho_air = rho_air,
        each mu_air = mu_air,
        each Pvacuum=Pvacuum,
        each GAMMA = GAMMA,
        each DELTA = DELTA,
        each BB = BB,
        each k_st = k_st,
        T_g_start_in=T_g_start_inlet_collector,
        T_g_start_out=T_g_start_inlet_collector,
        T_t_start_in=T_t_start_inlet_collector,
        T_t_start_out=T_t_start_outlet_collector,
        each Dext_g=Dext_g,
        each th_g=th_g,
        each Dext_t=Dext_t,
        each th_t=th_t) annotation (Placement(transformation(extent={{-30,-16},{14,34}})));

      ThermoCycle.Components.FluidFlow.Pipes.Flow1DimInc[Ns] flow1DimInc(
        redeclare each package Medium = Medium1,
        redeclare each final model Flow1DimIncHeatTransferModel =
            FluidHeatTransferModel,
        each N=N,
        each Nt = Nt,
        each A=A_lateral,
        each V=V_tube_int,
        each Mdotnom=Mdotnom,
        each Unom=Unom,
        each pstart=pstart,
        Tstart_inlet=Tstart_inlet_collector,
        Tstart_outlet=Tstart_outlet_collector,
        each steadystate=steadystate_T_fl,
        each Discretization=Discretization) annotation (Placement(transformation(
            extent={{-27.5,-31.5},{27.5,31.5}},
            rotation=90,
            origin={34.5,7.5})));

      ThermoCycle.Interfaces.Fluid.FlangeA InFlow(redeclare package Medium = Medium1) annotation (Placement(
            transformation(extent={{28,-110},{48,-90}}), iconTransformation(extent={{28,-110},{48,-90}})));

      ThermoCycle.Interfaces.Fluid.FlangeB OutFlow(redeclare package Medium = Medium1) annotation (Placement(
            transformation(extent={{30,106},{50,126}}), iconTransformation(extent={{30,106},{50,126}})));

      Modelica.Blocks.Interfaces.RealInput v_wind
        annotation (Placement(transformation(extent={{-84,74},{-44,114}}),
            iconTransformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-64,94})));
      Modelica.Blocks.Interfaces.RealInput Theta
        annotation (Placement(transformation(extent={{-96,4},{-56,44}}),
            iconTransformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-62,32})));
      Modelica.Blocks.Interfaces.RealInput Tamb
        annotation (Placement(transformation(extent={{-94,-26},{-54,14}}),
            iconTransformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-64,-20})));
      Modelica.Blocks.Interfaces.RealInput DNI
        annotation (Placement(transformation(extent={{-96,-56},{-56,-16}}),
            iconTransformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-64,-72})));

       /**************************************************** SUMMARY ********************************************/
    public
      record SummaryBase
         replaceable Arrays T_profile;
         record Arrays
          parameter Integer n;
          parameter Integer Ns;
          Modelica.SIunits.Temperature[Ns,n] T_fluid;
          Modelica.SIunits.Temperature[Ns,n] T_int_t;
          Modelica.SIunits.Temperature[Ns,n] T_t;
          Modelica.SIunits.Temperature[Ns,n] T_ext_t;
          Modelica.SIunits.Temperature[Ns,n] T_int_g;
          Modelica.SIunits.Temperature[Ns,n] T_g;
          Modelica.SIunits.Temperature[Ns,n] T_ext_g;
         end Arrays;
            Real Eta_solarCollector "Total efficiency of solar collector";
            Modelica.SIunits.HeatFlux Philoss
          "Heat Flux lost to the environment";
             Modelica.SIunits.Power Q_htf
          "Total heat through the termal heat transfer fluid flowing in the solar collector";
      end SummaryBase;
       replaceable record SummaryClass = SummaryBase;
       SummaryClass Summary( T_profile( n=N,Ns=Ns, T_fluid = T_fluid_,   T_int_t=T_int_t_,  T_t=T_t_, T_ext_t=T_ext_t_,  T_int_g=T_int_g_,  T_g=T_g_, T_ext_g=T_ext_g_), Eta_solarCollector=Eta_solarCollector_,Philoss=Philoss_,Q_htf=Q_htf_);
    protected
    Modelica.SIunits.Temperature T_fluid_[Ns,N];
    Modelica.SIunits.Temperature T_int_t_[Ns,N];
    Modelica.SIunits.Temperature T_t_[Ns,N];
    Modelica.SIunits.Temperature T_ext_t_[Ns,N];
    Modelica.SIunits.Temperature T_int_g_[Ns,N];
    Modelica.SIunits.Temperature T_g_[Ns,N];
    Modelica.SIunits.Temperature T_ext_g_[Ns,N];
    Real Eta_solarCollector_;
    Modelica.SIunits.HeatFlux Philoss_;
    Modelica.SIunits.Power Q_htf_;

    /*************************************************************/
    /************************* MODELLING *************************/
    /*************************************************************/
    equation

      // Temperature profiles
      for i in 1:Ns loop //Loop for each collector in series
          for k in 1:N loop // Loop for each cells in a given collector

            /* temperature profile of the working fluid*/
            T_fluid_[i,k] = flow1DimInc[i].Cells[k].T;

            /* temperature profile of the metal tube */
            T_int_t_[i,k] = solAbs[i].T_int_t[k];
            T_t_[i,k] = solAbs[i].T_t[k];
            T_ext_t_[i,k] = solAbs[i].T_ext_t[k];

            /* temperature profile of the glass envelope */
            T_int_g_[i,k] = solAbs[i].T_int_g[k];
            T_g_[i,k] = solAbs[i].T_g[k];
            T_ext_g_[i,k] = solAbs[i].T_ext_g[k];

          end for;

      end for;

      // Collector performance
      Eta_solarCollector_ = sum(solAbs[:].Eta_TOT);
      Philoss_ = sum(solAbs[:].Phi_loss);
      Q_htf_ = sum(flow1DimInc[:].Q_tot) "Total power absorbed by the fluid";

      // Element connections
      for i in 1:Ns loop
          connect(solAbs[i].wall_int, flow1DimInc[i].Wall_int) annotation (Line(
          points={{11.8,9},{20.45,9},{20.45,7.5},{21.375,7.5}},
          color={255,0,0},
          smooth=Smooth.None));
          connect(DNI, solAbs[i].DNI) annotation (Line(
          points={{-76,-36},{-36,-36},{-36,-10},{-27.36,-10}},
          color={0,0,127},
          smooth=Smooth.None));
          connect(v_wind, solAbs[i].v_wind) annotation (Line(
          points={{-64,94},{-34,94},{-34,30},{-26.92,30}},
          color={0,0,127},
          smooth=Smooth.None));
          connect(Tamb, solAbs[i].Tamb) annotation (Line(
          points={{-74,-6},{-46,-6},{-46,2.5},{-26.92,2.5}},
          color={0,0,127},
          smooth=Smooth.None));
          connect(Theta, solAbs[i].Theta) annotation (Line(
          points={{-76,24},{-38,24},{-38,12},{-26.92,12},{-26.92,14.5}},
          color={0,0,127},
          smooth=Smooth.None));
      end for;

      for i in 1:Ns-1 loop
        connect(flow1DimInc[i].OutFlow,flow1DimInc[i+1].InFlow);
      end for;

      connect(InFlow, flow1DimInc[1].InFlow) annotation (Line(
          points={{38,-100},{36,-100},{36,-15.4167},{34.5,-15.4167}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(flow1DimInc[Ns].OutFlow, OutFlow) annotation (Line(
          points={{34.2375,30.4167},{34.2375,116},{40,116}},
          color={0,0,255},
          smooth=Smooth.None));
                                                                                                          annotation (Dialog(group="Heat transfer", tab="General"),
                  Diagram(coordinateSystem(extent={{-80,-100},{100,100}},
              preserveAspectRatio=false),
                          graphics), Icon(coordinateSystem(extent={{-80,-100},{100,
                100}},    preserveAspectRatio=false),
                                          graphics={
                                          Bitmap(extent={{-96,118},{126,-100}}, fileName=
                  "modelica://ThermoCycle/Resources/Images/Avatar_SF.jpg"),
                                              Text(
              extent={{-44,106},{40,78}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name"),
            Text(
              extent={{-56,-32},{-30,-40}},
              lineColor={0,0,0},
              textString="Tamb[K]"),
            Text(
              extent={{-56,-82},{-26,-92}},
              lineColor={0,0,0},
              textString="DNI"),
            Text(
              extent={{-64,22},{-14,14}},
              lineColor={0,0,0},
              textString="Theta[rad]"),
            Text(
              extent={{-62,82},{-18,72}},
              lineColor={0,0,0},
              textString="v_wind [m/s]")}),
         Documentation(info="<HTML>

<p><big>The <b>SolarField_Forristal_Inc</b> is the same model of the <a href=\"modelica://ThermoCycle.Components.Units.Solar.SolarField_Forristal\">SolarField_Forristal</a> except that the heat transfer fluid flowing
in the collectors is modeled as an incompressible fluids.
</HTML>"));
    end SolarField_Forristal_Inc;
  end Component;
  annotation (uses(Modelica(version="3.2.1")));
end package_PSA_SFERAII;
