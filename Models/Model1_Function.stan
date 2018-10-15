functions{
  //ODE
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    //PARAMETERS
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    real k_IPTG = p[1];
    real kml0 = p[2];
    real kml = p[3];
    real theta_T = p[4];
    real theta_aTc = p[5];
    real n_aTc = p[6];
    real n_T = p[7];
    //real g_ml = p[8];
    real kmt0 = p[8];
    real kmt = p[9];
    real theta_L = p[10];
    real theta_IPTG = p[11];
    real n_IPTG = p[12];
    real n_L = p[13];
    //real g_mt = p[15];
    real k_pl = p[14];
    //real g_pl = p[17];
    real k_pt = p[15];
    //real g_pt = p[19];
    //Equations
    real dInd_dt[5];
    dInd_dt[1] = k_IPTG*(x_r[1]-y[1]);
    dInd_dt[2] = kml0+kml/(1+(y[5]/theta_T*1/(1+(x_r[2]/theta_aTc)^n_aTc))^n_T)-0.1386*y[2];
    dInd_dt[3] = kmt0+kmt/(1+(y[4]/theta_L*(1/(1+(y[1]/theta_IPTG)^n_IPTG))^n_L))-0.1386*y[3];
    dInd_dt[4] = k_pl*y[2]-0.0165*y[4];
    dInd_dt[5] = k_pt*y[3]-0.0165*y[5];
    //RESULTS
    return dInd_dt;
  }
  // When using expose_stan_functions() it allows to obtain a matrix with the results of the ODE over a time series
  matrix solve_coupled_ode(real[] ts, real[] y0, real[] p, real[] x_r){
    int x_i[0];
    return(to_matrix(integrate_ode_rk45(Toogle_one, y0,0,ts,p,x_r, x_i)));
  }
}
