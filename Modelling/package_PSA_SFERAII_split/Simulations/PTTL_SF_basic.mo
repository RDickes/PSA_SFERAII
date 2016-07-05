within package_PSA_SFERAII_split.Simulations;
model PTTL_SF_basic "SF + Mdot_source + SinkP"

  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot Supply(
    Mdot_0=5,
    redeclare package Medium = package_PSA_SFERAII_split.Media.Sytherm800,
    T_0=423.15)
    annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP Sink(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800, p0=1200000)
    annotation (Placement(transformation(extent={{78,44},{98,64}})));
  Component.SolarField_Forristal_Inc EuroTrough(
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
    Discretization=ThermoCycle.Functions.Enumerations.Discretizations.upwind_AllowFlowReversal,
    N=20,
    Unom=2000,
    redeclare model FluidHeatTransferModel =
        ThermoCycle.Components.HeatFlow.HeatTransfer.Constant,
    T_g_start_in=373.15,
    T_g_start_out=373.15,
    T_t_start_in=373.15,
    T_t_start_out=373.15,
    Tstart_inlet=373.15,
    Tstart_outlet=373.15,
    pstart=1000000)
    annotation (Placement(transformation(extent={{-12,-28},{44,36}})));

 Sources.FullDay_20160629.DNI DNI
    annotation (Placement(transformation(extent={{-90,-12},{-70,8}})));
  Sources.FullDay_20160629.theta theta
    annotation (Placement(transformation(extent={{-90,44},{-70,64}})));
  Sources.FullDay_20160629.V_wind v_wind
    annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
  Sources.FullDay_20160629.T_amb T_amb
    annotation (Placement(transformation(extent={{-90,16},{-70,36}})));
  Sources.FullDay_20160629.M_dot_htf m_dot_htf
    annotation (Placement(transformation(extent={{-90,-96},{-70,-76}})));
  Sources.FullDay_20160629.T_htf_su T_htf_su
    annotation (Placement(transformation(extent={{-90,-68},{-70,-48}})));
  ThermoCycle.Components.FluidFlow.Sensors.SensTp SensTex(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800)
    annotation (Placement(transformation(extent={{24,58},{38,72}})));
  ThermoCycle.Components.FluidFlow.Sensors.SensTp SensTsu(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800)
    annotation (Placement(transformation(extent={{-6,-92},{14,-72}})));
  Modelica.Blocks.Math.RealToInteger realToInteger
    annotation (Placement(transformation(extent={{-48,-40},{-36,-28}})));
  Sources.FullDay_20160629.Focus focus
    annotation (Placement(transformation(extent={{-90,-38},{-70,-18}})));
  Sources.FullDay_20160629.T_htf_ex t_htf_ex
    annotation (Placement(transformation(extent={{-8,74},{8,90}})));
  Component.Delta_T DeltaT
    annotation (Placement(transformation(extent={{52,78},{72,98}})));
equation
  connect(DNI.y, EuroTrough.DNI) annotation (Line(
      points={{-69,-2},{-54,-2},{-54,-10.72},{-8.26667,-10.72}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(theta.y, EuroTrough.Theta) annotation (Line(
      points={{-69,54},{-40,54},{-40,16.16},{-8.26667,16.16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(v_wind.y, EuroTrough.v_wind) annotation (Line(
      points={{-69,80},{-26,80},{-26,28},{-8.57778,28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(m_dot_htf.y, Supply.in_Mdot) annotation (Line(
      points={{-69,-86},{-34,-86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_amb.y, EuroTrough.Tamb) annotation (Line(
      points={{-69,26},{-44,26},{-44,2.72},{-8.26667,2.72}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_htf_su.y, Supply.in_T) annotation (Line(
      points={{-69,-58},{-28,-58},{-28,-86},{-28.2,-86}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(EuroTrough.OutFlow, SensTex.InFlow) annotation (Line(
      points={{26.5778,41.12},{24.6667,41.12},{24.6667,58.42},{31,58.42}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(Sink.flangeB, SensTex.InFlow) annotation (Line(
      points={{79.6,54},{58,54},{58,58.42},{31,58.42}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(Supply.flangeB,SensTsu. InFlow) annotation (Line(
      points={{-19,-92},{-8,-92},{-8,-91.4},{4,-91.4}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(SensTsu.InFlow, EuroTrough.InFlow) annotation (Line(
      points={{4,-91.4},{4,-91.7},{25.9556,-91.7},{25.9556,-30.56}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(realToInteger.y, EuroTrough.Focus) annotation (Line(
      points={{-35.4,-34},{-24,-34},{-24,-23.52},{-8.26667,-23.52}},
      color={255,127,0},
      smooth=Smooth.None));
  connect(focus.y, realToInteger.u) annotation (Line(
      points={{-69,-28},{-60,-28},{-60,-34},{-49.2,-34}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(t_htf_ex.y, DeltaT.EXP) annotation (Line(
      points={{8.8,82},{30,82},{30,83},{51.4,83}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(SensTex.T, DeltaT.SIM) annotation (Line(
      points={{36.6,69.2},{42,69.2},{42,93},{51.8,93}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),      graphics));
end PTTL_SF_basic;
