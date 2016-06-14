within ;
package package_PSA_SFERAII
  "package including all the models required for the PSA_SFERAII_project"
  annotation (uses(Modelica(version="3.2.1")));
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
  end Component;
end package_PSA_SFERAII;
