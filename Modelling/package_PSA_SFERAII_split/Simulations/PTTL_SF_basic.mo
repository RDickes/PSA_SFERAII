within package_PSA_SFERAII_split.Simulations;
model PTTL_SF_basic "SF + Mdot_source + SinkP"
  Component.SolarField_Forristal_Inc PTC_field(
    redeclare package Medium1 = package_PSA_SFERAII.Media.Sytherm800,
    eps1=1,
    eps2=1,
    eps3=1,
    eps4=1,
    eps5=1,
    eps6=1,
    rho_cl=0.931,
    Tau_g=0.92,
    Alpha_g=0.02,
    Eps_g=0.86,
    Alpha_t=0.7986,
    a1_IAM=4.11e-3,
    a2_IAM=5.513e-5,
    N=10,
    Ns=1,
    Nt=1,
    L=6*11.98,
    A_P=5.76,
    Dext_g=0.12,
    th_g=0.0025,
    rho_g=2230,
    Cp_g=900,
    lambda_g=1.14,
    Dext_t=0.07,
    th_t=0.002,
    rho_t=7990,
    Cp_t=515,
    lambda_t=18,
    Unom=100,
    redeclare model FluidHeatTransferModel =
        ThermoCycle.Components.HeatFlow.HeatTransfer.Ideal)
    "Field of Eurothrough  PTC"
    annotation (Placement(transformation(extent={{-24,-38},{36,22}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot Supply(Mdot_0=1,
      redeclare package Medium = package_PSA_SFERAII.Media.Sytherm800)
    annotation (Placement(transformation(extent={{-62,-86},{-42,-66}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP Sink(redeclare package
      Medium = package_PSA_SFERAII.Media.Sytherm800, p0=200000)
    annotation (Placement(transformation(extent={{70,62},{90,82}})));
equation
  connect(Supply.flangeB, PTC_field.InFlow) annotation (Line(
      points={{-43,-76},{15.3333,-76},{15.3333,-38}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(PTC_field.OutFlow, Sink.flangeB) annotation (Line(
      points={{16,26.8},{16,72},{71.6,72}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
            -100,-100},{100,100}}), graphics));
end PTTL_SF_basic;
