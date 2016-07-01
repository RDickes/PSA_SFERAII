within package_PSA_SFERAII_split.Simulations.BaseTests;
model test_syltherm800
  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot(redeclare
      package Medium = package_PSA_SFERAII_split.Media.Sytherm800, Mdot_0=1)
    annotation (Placement(transformation(extent={{-74,10},{-54,30}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP(redeclare package
      Medium = package_PSA_SFERAII_split.Media.Sytherm800)
    annotation (Placement(transformation(extent={{50,10},{70,30}})));
  ThermoCycle.Components.FluidFlow.Pipes.Flow1DimInc flow1DimInc(
    redeclare package Medium = package_PSA_SFERAII_split.Media.Sytherm800,
    Unom=10,
    pstart=100000,
    Tstart_inlet=283.15,
    Tstart_outlet=283.15)
    annotation (Placement(transformation(extent={{-14,8},{10,32}})));
equation
  connect(sourceMdot.flangeB, flow1DimInc.InFlow) annotation (Line(
      points={{-55,20},{-12,20}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(flow1DimInc.OutFlow, sinkP.flangeB) annotation (Line(
      points={{8,20.1},{30,20.1},{30,20},{51.6,20}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end test_syltherm800;
