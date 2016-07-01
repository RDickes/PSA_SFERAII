within package_PSA_SFERAII_split;
package UtilitiesSFERAII
  package Spline

    function MakeSpline
      extends Modelica.Icons.Function;
      input Integer n;
      input Real xp[n] "abscissa";
      input Real yp[n] "ordinates";
      input Real dydxl=Modelica.Constants.inf "left slope";
      input Real dydxr=Modelica.Constants.inf "right slope";
      output Real y2[n] "2nd derivatives";
    protected
      Real uh[n] "auxiliaey array";
      Real sig;
      Real p;
    algorithm
      if dydxl == Modelica.Constants.inf then
    //natural spline i.e. left curvature = 0
        y2[1]:= 0;
        uh[1]:= 0;
      else
    //given slope at the left border
        y2[1]:= -0.5;
        uh[1]:= +3/(xp[2] - xp[1])*((yp[2] - yp[1])/(xp[2] - xp[1]) - dydxl);
      end if;
    //decomposition of tridiagonal algorithm
      for i in 2:n - 1 loop
        sig:= (xp[i] - xp[i-1])/(xp[i+1] - xp[i-1]);
        p:= sig*y2[i-1] + 2;
        y2[i]:= (sig - 1)/p;
        uh[i]:= (6*((yp[i+1] - yp[i])/(xp[i+1] - xp[i]) -
                    (yp[i] - yp[i-1])/(xp[i] - xp[i-1]))/
                    (xp[i+1] - xp[i-1]) - sig*uh[i-1])/p;
      end for;
      if dydxr == Modelica.Constants.inf then
    //natural spline i.e. right curvature = 0
        y2[n]:= 0;
        uh[n]:= 0;
      else
    //given slope at the right border
        y2[n]:= +0.5;
        uh[n]:= -3/(xp[n] - xp[n-1])*((yp[n] - yp[n-1])/(xp[n] - xp[n-1]) - dydxr);
      end if;
    //backsubstitution of tridiagonal algorithm
      y2[n]:= (uh[n] - y2[n]*uh[n-1])/(y2[n]*y2[n-1] + 1);
      for i in 1:n - 1 loop
        y2[n-i]:= y2[n-i]*y2[n-i+1] + uh[n-i];
      end for;
    end MakeSpline;

    function EvalSpline
      extends Modelica.Icons.Function;
      input Real x;
      input Real xp[:] "abscissa";
      input Real yp[:] "ordinates";
      input Real y2[:] "2nd derivatives";
      output Real y;
    protected
      Integer n = size(xp,1);
      Integer i;
      Real h;
      Real a;
      Real b;
    algorithm
      i:= Bisection(x, xp);
      if i == 0 then
    //extrapolation below first point
        h:= (xp[2] - xp[1]);
        h:= (yp[2] - yp[1])/h - (2*y2[1] - (-1)*y2[2])*h/6;
        y:= yp[1] + h*(x - xp[1]);
      elseif i == n then
    //extrapolation above last point
        h:= (xp[n] - xp[n-1]);
        h:= (yp[n] - yp[n-1])/h - ((-1)*y2[n-1] - 2*y2[n])*h/6;
        y:= yp[n] + h*(x - xp[n]);
      else
    //evaluation of cubic spline
        h:= (xp[i+1] - xp[i]);
        a:= (xp[i+1] - x)/h;
        b:= (x - xp[i])/h;
        y:= a*yp[i] + b*yp[i+1] + ((a^3 - a)*y2[i] + (b^3 - b)*y2[i+1])*h^2/6;
      end if;
      annotation(derivative(noDerivative=xp, noDerivative=yp, noDerivative=y2)=Der1Spline);
    end EvalSpline;

    function Der1Spline
      extends Modelica.Icons.Function;
      input Real x;
      input Real xp[:] "abscissa";
      input Real yp[:] "ordinates";
      input Real y2[:] "2nd derivatives";
      input Real der_x;
      output Real der_y;
    protected
      Integer n = size(xp,1);
      Integer i;
      Real h;
      Real a;
      Real b;
    algorithm
      i:= Bisection(x, xp);
      if i == 0 then
    //extrapolation below first point
        h:= (xp[2] - xp[1]);
        der_y:= ((yp[2] - yp[1])/h - (2*y2[1] - (-1)*y2[2])*h/6)*der_x;
      elseif i == n then
    //extrapolation above last point
        h:= (xp[n] - xp[n-1]);
        der_y:= ((yp[n] - yp[n-1])/h - ((-1)*y2[n-1] - 2*y2[n])*h/6)*der_x;
      else
    //1st derivative of cubic spline
        h:= (xp[i+1] - xp[i]);
        a:= (xp[i+1] - x)/h;
        b:= (x - xp[i])/h;
        der_y:= ((yp[i+1] - yp[i])/h - ((3*a^2 - 1)*y2[i] - (3*b^2-1)*y2[i+1])*h/6)*der_x;
      end if;
    end Der1Spline;

    function Bisection
      extends Modelica.Icons.Function;
      input Real x;
      input Real xp[:] "abscissa";
      output Integer i;
    protected
      Integer n = size(xp,1);
      Integer ilo;
      Integer ihi;
    algorithm
      if x <= xp[1] then
    //below first point
        i:= 0;
      elseif x >= xp[n] then
    //above last point
        i:= n;
      else
    //bisectional search for the interval
        ilo:= 1;
        ihi:= n;
        while (ihi - ilo) > 1 loop
          i:= integer((ihi + ilo)/2);
          if (xp[i] > x) then
            ihi:= i;
          else
            ilo:= i;
          end if;
        end while;
    //xp[i] <= x < xp[i+1]
        i:= ilo;
      end if;
    end Bisection;

    model SplineFunction
      extends Modelica.Blocks.Interfaces.SISO;
      parameter Real table[:,2] "u data on first column, y data on y column";
    protected
      parameter Integer n = size(table,1) "Number of data points";
      parameter Real y2[n]= MakeSpline(n, table[:,1], table[:,2])
        "Spline derivative data";
    equation
      y = EvalSpline(u, table[:,1], table[:,2], y2);
      annotation (Icon(graphics={
            Line(points={{-78,82},{-78,-82}}, color={0,0,0}),
            Line(points={{-86,-72},{84,-72}}, color={0,0,0}),
            Line(points={{70,-64},{84,-72},{70,-78}}, color={0,0,0}),
            Line(points={{-84,68},{-78,82},{-72,68}}, color={0,0,0}),
            Line(points={{-70,-22},{-60,-14},{-46,-6},{-32,-2},{-18,-4},{-8,-10},{
                  2,-16},{14,-22},{28,-26},{38,-26},{48,-20},{56,-10},{62,2},{66,14}},
                color={0,0,0})}));
    end SplineFunction;

    model Test

      SplineFunction SplineFunction1(table=[0,300; 0.5,350; 0.8,400; 1,420])
        annotation (Placement(transformation(extent={{-40,20},{-20,40}}, rotation=
               0)));
      Modelica.Blocks.Sources.Ramp Ramp1(
        height=1,
        duration=0.9,
        offset=0,
        startTime=0.1)
                  annotation (Placement(transformation(extent={{-80,22},{-60,42}},
              rotation=0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
        PrescribedTemperature1
        annotation (Placement(transformation(extent={{20,0},{40,20}}, rotation=0)));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HeatCapacitor1(
                                                                 C=1)
        annotation (Placement(transformation(extent={{52,10},{72,30}}, rotation=0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
        PrescribedTemperature2
        annotation (Placement(transformation(extent={{20,-40},{40,-20}}, rotation=
               0)));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor HeatCapacitor2(
                                                                 C=1)
        annotation (Placement(transformation(extent={{52,-30},{72,-10}}, rotation=
               0)));
      SplineFunction SplineFunction2(table=[0,300; 1,350])
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}},
              rotation=0)));
    equation
      connect(Ramp1.y, SplineFunction1.u)
        annotation (Line(points={{-59,32},{-42,32},{-42,30}}, color={0,0,127}));
      connect(SplineFunction1.y, PrescribedTemperature1.T)
        annotation (Line(points={{-19,30},{8,30},{8,10},{18,10}}, color={0,0,127}));
      connect(PrescribedTemperature1.port, HeatCapacitor1.port)
        annotation (Line(points={{40,10},{62,10}}, color={191,0,0}));
      connect(PrescribedTemperature2.port,HeatCapacitor2. port)
        annotation (Line(points={{40,-30},{62,-30}}, color={191,0,0}));
      connect(SplineFunction2.y, PrescribedTemperature2.T)
        annotation (Line(points={{-19,-30},{18,-30}}, color={0,0,127}));
      connect(Ramp1.y, SplineFunction2.u) annotation (Line(points={{-59,32},{-54,32},
              {-54,-30},{-42,-30}}, color={0,0,127}));
      annotation (Diagram(graphics));
    end Test;
    annotation (uses(Modelica(version="3.1")),
      preferedView="info",
      Documentation(info="<html>
<p>
This package implements cubic spline interpolation as described in<br>
Numerical Recipes in FORTRAN 77, Chapter 3.3 Cubic Spline Interpolation
</p>
<p>
Copyright &copy; 2005 <a href=\"http://www.haumer.at\">Anton Haumer</a>.
</p>
</html>"));
  end Spline;
  annotation ();
end UtilitiesSFERAII;
