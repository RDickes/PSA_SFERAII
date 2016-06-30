within package_PSA_SFERAII_split.Simulations;
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
    annotation (Placement(transformation(extent={{-32,-26},{22,38}})));

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
      points={{-75.3,73},{-55.65,73},{-55.65,36.08},{-27.2,36.08}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(theta_input.y, solarField_Forristal_Inc.Theta) annotation (Line(
      points={{-75.2,42},{-62,42},{-62,16.24},{-26.6,16.24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_amb_input.y, solarField_Forristal_Inc.Tamb) annotation (Line(
      points={{-77.3,9},{-57.65,9},{-57.65,-0.4},{-27.2,-0.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(DNI_input.y, solarField_Forristal_Inc.DNI) annotation (Line(
      points={{-77.3,-21},{-57.65,-21},{-57.65,-17.04},{-27.2,-17.04}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceMdot.flangeB, solarField_Forristal_Inc.InFlow) annotation (
      Line(
      points={{-63,-72},{2,-72},{2,-26},{3.4,-26}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(solarField_Forristal_Inc.OutFlow, sinkP.flangeB) annotation (Line(
      points={{4,43.12},{4,80},{61.6,80}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
            -100,-100},{100,100}}), graphics));
end Basic1;
