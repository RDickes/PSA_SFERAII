within package_PSA_SFERAII_split.Component;
model SolarField_Forristal_Inc
  "PSA_SFERAII - Solar field model with collectors based on Forristal model for incompressible fluids"

/*************************************************************/
/************************* 1. FLUID **************************/
/*************************************************************/
replaceable package Medium1 = ThermoCycle.Media.DummyFluid constrainedby
    Modelica.Media.Interfaces.PartialMedium                                                                      annotation (choicesAllMatching = true);

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
parameter Modelica.SIunits.Length Dext_g = 0.12 "External glass diameter [m] "
                                                                              annotation (Dialog(group="Properties of the glass envelope", tab="General"));
parameter Modelica.SIunits.Length th_g = 0.0025 "Glass thickness [m] "
                                                                      annotation (Dialog(group="Properties of the glass envelope", tab="General"));
parameter Boolean GlassUD = true
    "if true, constant properties defined by the user"                               annotation (Dialog(group="Properties of the glass envelope", tab="General"));
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
parameter Modelica.SIunits.Length Dext_t =  0.07 "External diameter tube [m] "
                                                                              annotation (Dialog(group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.Length th_t =  0.002 "tube thickness [m] "
                                                                     annotation (Dialog(group="Properties of the tube receiver", tab="General"));
final parameter Modelica.SIunits.Length D_int_t= Dext_t - 2*th_t
    "internal diameter of the metal tube [m] ";
parameter Boolean TubeUD = true
    "if true, constant properties defined by the user"                              annotation (Dialog(group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.Density rho_t = 8000 "tube density [kg/m3]"  annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.SpecificHeatCapacity Cp_t = 500
    "Specific heat capacity of the tube [J/kg.K] "                                                         annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
parameter Modelica.SIunits.ThermalConductivity lambda_t = 54
    "Thermal conductivity of the tube [W/m.K] "                                                          annotation (Dialog(enable=TubeUD,group="Properties of the tube receiver", tab="General"));
replaceable parameter Material.TubeReceiver.MaterialBase TubeMaterial constrainedby
    package_PSA_SFERAII_split.Component.Material.TubeReceiver.MaterialBase
    "Select tube material for the receiver if not user defined"                                     annotation (choicesAllMatching=true, Dialog(enable=not
                                                                                                (TubeUD),group="Properties of the tube receiver", tab="General"));

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
parameter Real GAMMA= 1.39 "Ratio of specific heats for the annulus gas [-]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
parameter Real DELTA = 3.53e-8 "Molecular diameter for the annulus gas [cm]" annotation (Dialog(group="Vacuum properties: between metal and glass envelope", tab="General"));
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
parameter Modelica.SIunits.MassFlowRate Mdotnom= 0.5 "Total nominal Mass flow";
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

    package_PSA_SFERAII_split.Component.SolAbsForristal[Ns] solAbs(
    each eps1=eps1,
    each eps2=eps2,
    each eps3=eps3,
    each rho_cl=rho_cl,
    each eps4=eps4,
    each eps5=eps5,
    each eps6=eps6,
    each Tau_g=Tau_g,
    each Alpha_g=Alpha_g,
    each Eps_g=Eps_g,
    each Alpha_t=Alpha_t,
    each a1_IAM=a1_IAM,
    each a2_IAM=a2_IAM,
    each N=N,
    each L=L,
    each A_P=A_P,
    each GlassUD=GlassUD,
    each rho_g_ud=rho_g_ud,
    each Cp_g_ud=Cp_g_ud,
    each lambda_g_ud=lambda_g_ud,
    each GlassMaterial=GlassMaterial,
    each TubeUD=TubeUD,
    each rho_t_ud=rho_t_ud,
    each Cp_t_ud=Cp_t_ud,
    each lambda_t_ud=lambda_t_ud,
    each TubeMaterial=TubeMaterial,
    each Patm=Patm,
    each Pr=Pr,
    each k_air=k_air,
    each rho_air=rho_air,
    each mu_air=mu_air,
    each Pvacuum=Pvacuum,
    each GAMMA=GAMMA,
    each DELTA=DELTA,
    each BB=BB,
    each k_st=k_st,
    T_g_start_in=T_g_start_inlet_collector,
    T_g_start_out=T_g_start_inlet_collector,
    T_t_start_in=T_t_start_inlet_collector,
    T_t_start_out=T_t_start_outlet_collector,
    each Dext_g=Dext_g,
    each th_g=th_g,
    each Dext_t=Dext_t,
    each th_t=th_t)
    annotation (Placement(transformation(extent={{-30,-16},{14,34}})));

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
        Modelica.SIunits.HeatFlux Philoss "Heat Flux lost to the environment";
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
