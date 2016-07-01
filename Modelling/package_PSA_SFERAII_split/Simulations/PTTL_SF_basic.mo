within package_PSA_SFERAII_split.Simulations;
model PTTL_SF_basic "SF + Mdot_source + SinkP"


  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot Supply(
    Mdot_0=5,
    redeclare package Medium = package_PSA_SFERAII_split.Media.Sytherm800,
    T_0=423.15)
    annotation (Placement(transformation(extent={{-38,-94},{-18,-74}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP Sink(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800, p0=1200000)
    annotation (Placement(transformation(extent={{70,52},{90,72}})));
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
    T_g_start_in=373.15,
    T_g_start_out=373.15,
    T_t_start_in=373.15,
    T_t_start_out=373.15,
    Unom=2000,
    Tstart_inlet=373.15,
    Tstart_outlet=373.15,
    pstart=1000000,
    redeclare model FluidHeatTransferModel =
        ThermoCycle.Components.HeatFlow.HeatTransfer.Constant)
    annotation (Placement(transformation(extent={{-14,-22},{42,42}})));
  ExpData.FullDay_20160629.DNI DNI
    annotation (Placement(transformation(extent={{-90,-18},{-70,2}})));
  ExpData.FullDay_20160629.theta theta
    annotation (Placement(transformation(extent={{-90,44},{-70,64}})));
  ExpData.FullDay_20160629.V_wind v_wind
    annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
  ExpData.FullDay_20160629.T_amb T_amb
    annotation (Placement(transformation(extent={{-90,16},{-70,36}})));
  ExpData.FullDay_20160629.M_dot_htf m_dot_htf
    annotation (Placement(transformation(extent={{-90,-84},{-70,-64}})));
  ExpData.FullDay_20160629.T_htf_su T_htf_su
    annotation (Placement(transformation(extent={{-90,-56},{-70,-36}})));
equation
  connect(Supply.flangeB, EuroTrough.InFlow) annotation (Line(
      points={{-19,-84},{22,-84},{22,-22},{22.7111,-22}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(DNI.y, EuroTrough.DNI) annotation (Line(
      points={{-69,-8},{-42,-8},{-42,-13.04},{-9.02222,-13.04}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(theta.y, EuroTrough.Theta) annotation (Line(
      points={{-69,54},{-40,54},{-40,20.24},{-8.4,20.24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(v_wind.y, EuroTrough.v_wind) annotation (Line(
      points={{-69,80},{-26,80},{-26,40.08},{-9.02222,40.08}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(m_dot_htf.y, Supply.in_Mdot) annotation (Line(
      points={{-69,-74},{-34,-74},{-34,-78}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_amb.y, EuroTrough.Tamb) annotation (Line(
      points={{-69,26},{-44,26},{-44,3.6},{-9.02222,3.6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_htf_su.y, Supply.in_T) annotation (Line(
      points={{-69,-46},{-28,-46},{-28,-78},{-28.2,-78}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(EuroTrough.OutFlow, Sink.flangeB) annotation (Line(
      points={{23.3333,47.12},{23.6667,47.12},{23.6667,62},{71.6,62}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),      graphics));
end PTTL_SF_basic;
