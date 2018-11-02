functions{
  //ODE
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
   // Inputs
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    // Parameters
    real k_in_IPTG = p[1];
    real k_out_IPTG = p[2];
    real k_in_aTc = p[3];
    real k_out_aTc = p[4];
    real k_L_pm0 = p[5];
    real k_L_pm = p[6];
    real theta_T = p[7];
    real theta_aTc = p[8];
    real n_aTc = p[9];
    real n_T = p[10];
    real k_T_pm0 = p[11];
    real k_T_pm = p[12];
    real theta_L = p[13];
    real theta_IPTG = p[14];
    real n_IPTG = p[15];
    real n_L = p[16];
    
    //Equations
    real dInd_dt[4];
    
    if (x_r[1]>y[1]){
      dInd_dt[1]=k_in_IPTG*(x_r[1]-y[1]);
    }
    else{
      dInd_dt[1]=k_out_IPTG*(x_r[1]-y[1]);
    }
    
    if (x_r[2]>y[2]){
      dInd_dt[2]=k_in_aTc*(x_r[2]-y[2]);
    }
    else{
      dInd_dt[2]=k_out_aTc*(x_r[2]-y[2]);
    }
    
    dInd_dt[3] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+(y[4]/theta_T*1/(1+(y[2]/theta_aTc)^n_aTc))^n_T))))-0.0165*y[3];
    dInd_dt[4] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+(y[3]/theta_L*1/(1+(y[1]/theta_IPTG)^n_IPTG))^n_L))))-0.0165*y[4];
    // RESULTS
    return dInd_dt;
    
  }
  matrix solve_coupled_ode(real[] ts, real[] y0ss, real[] p, real[] iss, real[] inp, real[] x_r, int[] x_i){
    real steady[num_elements(ts),4];
    real simul[num_elements(ts),4];
    real y0[4];
    
    steady = integrate_ode_rk45(Toogle_one, y0ss,0,ts,p,iss, x_i, 1e-8, 1e-8, 1e7);
    
    y0 = steady[num_elements(steady[,1]),];
    
    simul = integrate_ode_rk45(Toogle_one, y0,0,ts,p,inp, x_i, 1e-8, 1e-8, 1e7);
    

    // print(simul[num_elements(simul[,1]),])
    return(to_matrix(simul));
    
  }
}
