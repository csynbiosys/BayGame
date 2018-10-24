functions{
  //ODE
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    //PARAMETERS
    
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    
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
  matrix solve_coupled_ode(real[] ts, real[] y0, real[] p, real[] x_r, int[] x_i, int[] sp, real[] inputs){
    int maxtime = num_elements(ts);
    int Nsp = num_elements(sp);
    int Nevents = num_elements(sp)-1;
    int Neq = 4;
    
    // matrix[maxtime,Neq] total;
    real final[maxtime,Neq];
    real initialV[Neq];
    int i;
    
    
    initialV = y0;
    i = 1;
    
    for (q in 1:Nevents){
      int itp = sp[q];  // General way to extract the initial time points of each event
      int lts = num_elements(ts[(sp[q]+1):sp[q+1]]);  // General way to define the number of elements in each event series
      real input = inputs[q]; // General way to extract the input values
      real Tevent[lts] = ts[(sp[q]+1):sp[q+1]];  // General way to extract the times of each event
      
      real part1[lts,Neq];
      real temp[lts,Neq];
      
      
      if (q == 1){part1 = integrate_ode_rk45(Toogle_one, initialV,itp,ts[(sp[q]+1):sp[q+1]],p,to_array_1d(inputs[i:(i+1)]), x_i);}
      else{part1 = integrate_ode_rk45(Toogle_one, initialV,(itp-1),ts[(sp[q]+1):sp[q+1]],p,to_array_1d(inputs[i:(i+1)]), x_i);}
      
      initialV = part1[lts];
      i=i+2;
      
      for (y in (itp+1):(itp+lts)){
        
        // total[(y),]=to_matrix(part1)[(y-itp),];
        final[(y),]=(part1)[(y-itp),];
      };
      
    };
    
    return(to_matrix(final));
  }
}
