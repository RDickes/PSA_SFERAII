within package_PSA_SFERAII_split.Component.Material.GlassEnvelope;
record MaterialBase "SolidMaterial example for the glass envelope in PTCs"

parameter Real a0_rho_g = 2707 "rho_g = a0_rho_g + a1_rho_g*T + a2_rho_g*T^2";
parameter Real a1_rho_g = 0 "rho_g = a0_rho_g + a1_rho_g*T + a2_rho_g*T^2";
parameter Real a2_rho_g = 0 "rho_g = a0_rho_g + a1_rho_g*T + a2_rho_g*T^2";

parameter Real a0_cp_g = 2707 "cp_g = a0_cp_g + a1_cp_g*T + a2_cp_g*T^2";
parameter Real a1_cp_g = 0 "cp_g = a0_cp_g + a1_cp_g*T + a2_cp_g*T^2";
parameter Real a2_cp_g = 0 "cp_g = a0_cp_g + a1_cp_g*T + a2_cp_g*T^2";

parameter Real a0_lambda_g = 2707
    "lambda_g = a0_lambda_g + a1_lambda_g*T + a2_lambda_g*T^2";
parameter Real a1_lambda_g = 0
    "lambda_g = a0_lambda_g + a1_lambda_g*T + a2_lambda_g*T^2";
parameter Real a2_lambda_g = 0
    "lambda_g = a0_lambda_g + a1_lambda_g*T + a2_lambda_g*T^2";

end MaterialBase;
