functions{
  //ODE
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    //PARAMETERS
    real u_aTc = x_r[1];
    real u_IPTG = x_r[2];
    real k_in_aTc = p[1];
    real k_out_aTc = p[2];
    real k_in_IPTG = p[3];
    real k_out_IPTG = p[4];
    real kml0 = p[5];
    real kml = p[6];
    real theta_T = p[7];
    real theta_aTc = p[8];
    real n_aTc = p[9];
    real n_T = p[10];
    // real g_ml = p[11];
    real kmt0 = p[11];
    real kmt = p[12];
    real theta_L = p[13];
    real theta_IPTG = p[14];
    real n_IPTG = p[15];
    real n_L = p[16];
    // real g_mt = p[18];
    real k_pl = p[17];
    // real g_pl = p[20];
    real k_pt = p[18];
    // real g_pt = p[22];
    //Equations
    real dInd_dt[6];
    if (x_r[1]>y[1]){
      dInd_dt[1]=k_in_aTc*(x_r[1]-y[1]);
    }
    else{
      dInd_dt[1]=k_out_aTc*(x_r[1]-y[1]);
    }
    if (x_r[2]>y[2]){
      dInd_dt[2]=k_in_IPTG*(x_r[2]-y[2]);
    }
    else{
      dInd_dt[2]=k_out_IPTG*(x_r[2]-y[2]);
    }
    dInd_dt[3] = kml0+(kml/(1+(y[6]/theta_T*1/(1+(y[1]/theta_aTc)^n_aTc))^n_T))-0.1386*y[3];
    dInd_dt[4] = kmt0+(kmt/(1+(y[5]/theta_L*1/(1+(y[4]/theta_IPTG)^n_IPTG))^n_L))-0.1386*y[4];
    dInd_dt[5] = k_pl*y[3]-0.0165*y[5];
    dInd_dt[6] = k_pt*y[4]-0.0165*y[6];
    // RESULTS
    return dInd_dt;
  }
  matrix solve_coupled_ode(real[] ts, real[] y0, real[] p, real[] x_r){
    int x_i[0];
    return(to_matrix(integrate_ode_rk45(Toogle_one, y0, 0, ts, p, x_r, x_i)));
  }
}
