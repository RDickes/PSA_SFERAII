within package_PSA_SFERAII_split.Media;
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
