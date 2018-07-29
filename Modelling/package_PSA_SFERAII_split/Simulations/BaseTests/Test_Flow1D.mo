within package_PSA_SFERAII_split.Simulations.BaseTests;
model Test_Flow1D
  replaceable package Medium =
     package_PSA_SFERAII_split.Media.Sytherm800                                                  constrainedby
    Modelica.Media.Interfaces.PartialMedium;

  ThermoCycle.Components.FluidFlow.Pipes.Flow1DimInc flow1DimInc(redeclare
      package Medium =                                                                Medium,
    Nt=1,
    N=10,
    V=0.26,
    Mdotnom=9.5,
    A=15.35,
    Unom=1200,
    Discretization=ThermoCycle.Functions.Enumerations.Discretizations.upwind_AllowFlowReversal,
    pstart=400000,
    Tstart_inlet=473.15,
    Tstart_outlet=513.15,
    redeclare model Flow1DimIncHeatTransferModel =
        ThermoCycle.Components.HeatFlow.HeatTransfer.SinglePhase (redeclare
          model LiquidCorrelation =
            ThermoCycle.Components.HeatFlow.HeatTransfer.SinglePhaseCorrelations.DittusBoelter1930
            (d_h=0.066)))
    annotation (Placement(transformation(extent={{-34,-20},{16,30}})));

  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot( redeclare
      package Medium =                                                                          Medium,
    Mdot_0=9.5,
    T_0=473.15)
    annotation (Placement(transformation(extent={{-84,-2},{-64,18}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP( redeclare package
      Medium =                                                                        Medium)
    annotation (Placement(transformation(extent={{60,8},{80,28}})));
  ThermoCycle.Components.HeatFlow.Sources.HeatSource heatSource(N=10,
    HeatFlow=true,
    A=15.35)
    annotation (Placement(transformation(extent={{-30,30},{12,58}})));
  Modelica.Blocks.Sources.Constant const(k=7.5e5)
    annotation (Placement(transformation(extent={{-68,56},{-48,76}})));
equation
  connect(flow1DimInc.OutFlow, sinkP.flangeB) annotation (Line(
      points={{11.8333,5.20833},{36,5.20833},{36,18},{61.6,18}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(sourceMdot.flangeB, flow1DimInc.InFlow) annotation (Line(
      points={{-65,8},{-50,8},{-50,5},{-29.8333,5}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(heatSource.thermalPort, flow1DimInc.Wall_int) annotation (Line(
      points={{-9.21,38.26},{-9,38.26},{-9,15.4167}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(const.y, heatSource.Phi) annotation (Line(
      points={{-47,66},{-9,66},{-9,49.6}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end Test_Flow1D;
