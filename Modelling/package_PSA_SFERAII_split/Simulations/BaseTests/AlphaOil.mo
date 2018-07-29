within package_PSA_SFERAII_split.Simulations.BaseTests;
model AlphaOil
  import Modelica.SIunits;
  import Modelica.Constants;
  replaceable package Medium =
     package_PSA_SFERAII_split.Media.Sytherm800                                                  constrainedby
    Modelica.Media.Interfaces.PartialMedium;

parameter SIunits.Temperature T_su_oil=200+273.15;
parameter SIunits.Temperature T_ex_oil=240+273.15;
parameter SIunits.MassFlowRate m_oil = 9.5;
parameter SIunits.AbsolutePressure p_oil = 1e5;

/* Geometries parameter */
parameter SIunits.Length LL= 71.88;
parameter SIunits.Length dd_int = 0.068;
final parameter SIunits.Area A_l = LL * Constants.pi * dd_int;
final parameter SIunits.Area A_cross = dd_int^2/4 * Constants.pi;
final parameter SIunits.Volume V_tube = A_cross * LL;
/* Variables */
SIunits.Temperature T_avg;
SIunits.SpecificHeatCapacity cp_oil;
SIunits.Density rho_oil;
SIunits.DynamicViscosity mu_oil;
SIunits.ThermalConductivity k_oil;
SIunits.Velocity v_oil;
SIunits.Power   Q_oil;

Medium.ThermodynamicState  fluidState;

SIunits.ReynoldsNumber Re(min=0);
SIunits.PrandtlNumber Pr(min=0);

SIunits.CoefficientOfHeatTransfer alpha_oil;

/* Correlation DittusBolter 1930*/
parameter Real a_DB = 0.023;
parameter Real b_DB = 0.8;
parameter Real c_DB = 0.4;
SIunits.NusseltNumber Nu_DB;
SIunits.CoefficientOfHeatTransfer alpha_DB;

equation
  T_avg = (T_su_oil + T_ex_oil)/2;
  fluidState = Medium.setState_pT(p_oil, T_avg);

  cp_oil = Medium.specificHeatCapacityCp(fluidState);
  rho_oil = Medium.density(fluidState);
  mu_oil = Medium.dynamicViscosity(fluidState);
  k_oil = Medium.thermalConductivity(fluidState);

  Q_oil = m_oil * cp_oil * (T_ex_oil - T_su_oil);
  alpha_oil = Q_oil / ((T_ex_oil - T_su_oil)*A_l);

  v_oil = m_oil / (rho_oil * A_cross);
  Re = v_oil * rho_oil * dd_int / mu_oil;
  Pr = cp_oil * mu_oil /  k_oil;
  Nu_DB =  a_DB * Re^b_DB * Pr^c_DB;
  alpha_DB  = Nu_DB * k_oil / dd_int;

end AlphaOil;
