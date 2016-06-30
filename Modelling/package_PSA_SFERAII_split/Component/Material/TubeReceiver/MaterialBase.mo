within package_PSA_SFERAII_split.Component.Material.TubeReceiver;
record MaterialBase "Material example for the tube receiver in the PTCs"

parameter Real a0_rho_t = 2707 "rho_t = a0_rho_t + a1_rho_t*T + a2_rho_t*T^2";
parameter Real a1_rho_t = 0 "rho_t = a0_rho_t + a1_rho_t*T + a2_rho_t*T^2";
parameter Real a2_rho_t = 0 "rho_t = a0_rho_t + a1_rho_t*T + a2_rho_t*T^2";

parameter Real a0_cp_t = 2707 "cp_t = a0_cp_t + a1_cp_t*T + a2_cp_t*T^2";
parameter Real a1_cp_t = 0 "cp_t = a0_cp_t + a1_cp_t*T + a2_cp_t*T^2";
parameter Real a2_cp_t = 0 "cp_t = a0_cp_t + a1_cp_t*T + a2_cp_t*T^2";

parameter Real a0_lambda_t = 2707
    "lambda_t = a0_lambda_t + a1_lambda_t*T + a2_lambda_t*T^2";
parameter Real a1_lambda_t = 0
    "lambda_t = a0_lambda_t + a1_lambda_t*T + a2_lambda_t*T^2";
parameter Real a2_lambda_t = 0
    "lambda_t = a0_lambda_t + a1_lambda_t*T + a2_lambda_t*T^2";

end MaterialBase;
