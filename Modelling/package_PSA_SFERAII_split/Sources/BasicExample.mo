within package_PSA_SFERAII_split.Sources;
package BasicExample

  model DNI
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end DNI;

  model theta
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end theta;

  model V_wind
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end V_wind;

  model T_amb
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end T_amb;

  model T_htf_su
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end T_htf_su;

  model M_dot_htf
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end M_dot_htf;

  model P_htf
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end P_htf;

  model Focus
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end Focus;

  model T_htf_ex
    extends Spline.SplineBuilder(
        table=[1,1]);
  equation

    annotation (experiment(
        StartTime=-500,
        StopTime=35000,
        NumberOfIntervals=5000), experimentSetupOutput);
  end T_htf_ex;
end BasicExample;
