within package_PSA_SFERAII_split.Component;
model Delta_T

  Modelica.Blocks.Interfaces.RealInput SIM
    annotation (Placement(transformation(extent={{-122,30},{-82,70}})));
  Modelica.Blocks.Interfaces.RealInput EXP
    annotation (Placement(transformation(extent={{-126,-70},{-86,-30}})));
  Modelica.Blocks.Interfaces.RealOutput Delta
    annotation (Placement(transformation(extent={{92,-26},{146,28}})));

equation
  Delta = SIM-273.15 - EXP;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics));
end Delta_T;
