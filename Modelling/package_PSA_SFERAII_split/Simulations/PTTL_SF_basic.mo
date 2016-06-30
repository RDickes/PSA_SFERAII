within package_PSA_SFERAII_split.Simulations;
model PTTL_SF_basic "SF + Mdot_source + SinkP"


  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot Supply(
    Mdot_0=5,
    redeclare package Medium = package_PSA_SFERAII_split.Media.Sytherm800,
    T_0=423.15)
    annotation (Placement(transformation(extent={{-62,-86},{-42,-66}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP Sink(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800, p0=1200000)
    annotation (Placement(transformation(extent={{70,62},{90,82}})));
  Modelica.Blocks.Sources.Step v_wind_input(height=0)
    annotation (Placement(transformation(extent={{-80,76},{-66,90}})));
  Modelica.Blocks.Sources.Step theta_input(height=0)
    annotation (Placement(transformation(extent={{-82,18},{-66,34}})));
  Modelica.Blocks.Sources.Step T_amb_input(height=298)
    annotation (Placement(transformation(extent={{-82,-14},{-68,0}})));
  Modelica.Blocks.Sources.Step DNI_input(height=1000)
    annotation (Placement(transformation(extent={{-82,-44},{-68,-30}})));
  Component.SolarField_Forristal_Inc solarField_Forristal_Inc(
    redeclare package Medium1 = package_PSA_SFERAII_split.Media.Sytherm800,
    Ns=1,
    Nt=1,
    L=6*11.8,
    A_P=5.76,
    GlassUD=false,
    redeclare package_PSA_SFERAII_split.Component.Material.GlassEnvelope.Pyrex
      GlassMaterial,
    TubeUD=false,
    redeclare
      package_PSA_SFERAII_split.Component.Material.TubeReceiver.SLsteel_316
      TubeMaterial,
    T_g_start_in=373.15,
    T_g_start_out=373.15,
    T_t_start_in=373.15,
    T_t_start_out=373.15,
    Tstart_inlet=373.15,
    Tstart_outlet=373.15,
    pstart=1000000,
    Discretization=ThermoCycle.Functions.Enumerations.Discretizations.upwind_AllowFlowReversal)
    annotation (Placement(transformation(extent={{-14,-22},{42,42}})));
equation
  connect(v_wind_input.y, solarField_Forristal_Inc.v_wind) annotation (Line(
      points={{-65.3,83},{-39.65,83},{-39.65,40.08},{-9.02222,40.08}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(theta_input.y, solarField_Forristal_Inc.Theta) annotation (Line(
      points={{-65.2,26},{-38,26},{-38,20.24},{-8.4,20.24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_amb_input.y, solarField_Forristal_Inc.Tamb) annotation (Line(
      points={{-67.3,-7},{-39.65,-7},{-39.65,3.6},{-9.02222,3.6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(DNI_input.y, solarField_Forristal_Inc.DNI) annotation (Line(
      points={{-67.3,-37},{-40.65,-37},{-40.65,-13.04},{-9.02222,-13.04}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Supply.flangeB, solarField_Forristal_Inc.InFlow) annotation (Line(
      points={{-43,-76},{22,-76},{22,-22},{22.7111,-22}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(solarField_Forristal_Inc.OutFlow, Sink.flangeB) annotation (Line(
      points={{23.3333,47.12},{23.6667,47.12},{23.6667,72},{71.6,72}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),      graphics));
end PTTL_SF_basic;
