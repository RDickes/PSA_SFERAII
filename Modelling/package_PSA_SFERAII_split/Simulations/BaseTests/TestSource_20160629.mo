within package_PSA_SFERAII_split.Simulations.BaseTests;
model TestSource_20160629
  Sources.FullDay_20160629.DNI dNI
    annotation (Placement(transformation(extent={{-34,46},{-14,66}})));
  Sources.FullDay_20160629.T_htf_su t_htf_su
    annotation (Placement(transformation(extent={{-46,4},{-26,24}})));
  Sources.FullDay_20160629.M_dot_htf m_dot_htf
    annotation (Placement(transformation(extent={{-50,-46},{-30,-26}})));
  annotation (experiment(StopTime=55000), __Dymola_experimentSetupOutput);
end TestSource_20160629;
