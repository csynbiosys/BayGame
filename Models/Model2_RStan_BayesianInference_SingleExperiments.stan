// RStan model containing the second inducer exchange model for the Toggle Switch contained in the Lugagne paper

functions{
  // Function containing the ODE to be used for the inference
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
    // }

    // RESULTS
    return dInd_dt;
  }
  
  // Function type vector containing the equations where the root needs to be calculated for the steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    vector[2] alpha;
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
    
    // Equations
    alpha[1] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+((init[2]/theta_T)*1/(1+(x_r[2]/theta_aTc)^n_aTc))^n_T))))/0.0165;
    alpha[2] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+((init[1]/theta_L)*1/(1+(x_r[1]/theta_IPTG)^n_IPTG))^n_L))))/0.0165;
    // Results
    return alpha;
  }
  
}

data {
  // INPUTS
  int<lower = 0> tsl; // Length of time vector
  real ts[tsl]; // Total time vector
  int tsmax; // Last element of the time vector
  real preIPTG; // Pres values
  real preaTc;
  int Nsp; // Number of event switching points (including initial and final)
  int evnT[Nsp]; // Event switching points (including initial and final)
  real inputs[(Nsp-1)*2]; // Inputs as IPTG, aTc, IPTG, aTc, ...
  real IPTG[Nsp-1]; // Input values at each event
  real aTc[Nsp-1];
  
  // OBSERVABLES
  int<lower=0> stsl; // Number of sampling times
  int sts[stsl]; // Sampling times
  real GFPmean[stsl]; // Mean and standard deviations for proteins
  real RFPmean[stsl];
  real<lower=0> GFPstd[stsl];
  real<lower=0> RFPstd[stsl];
  
  int tonil;
  real toni[tonil];

}

transformed data {
  int nParms = 16; // Number of parameters of the model
  int Neq = 4; // Total number of equations of the model
  int Nevents = Nsp-1; // Number of events
  int x_i[0]; // Empty x_i object (needs to be deffined)
  real x_r[(Nevents)*2]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq-1]; // Initial experimental values for the calculation of the steady state ordered as IPTG, aTc
  real pre[2]; // Input values during the ON incubation ordered as IPTG, aTc
  
  ivss[1] = RFPmean[1];
  ivss[2] = GFPmean[1];
  
  pre[1] = preIPTG;
  pre[2] = preaTc;
  
}

parameters {
    // Parameters to be infered in the model
    real<lower=0> k_in_IPTG_raw;
    real<lower=0> k_out_IPTG_raw;
    real<lower=0> k_in_aTc_raw;
    real<lower=0> k_out_aTc_raw;
    real<lower=0> k_L_pm0_raw;
    real<lower=0> k_L_pm_raw;
    real<lower=0> theta_T_raw;
    real<lower=0,upper=100> theta_aTc_raw;
    real<lower=0> n_aTc_raw;
    real<lower=0> n_T_raw;
    real<lower=0> k_T_pm0_raw;
    real<lower=0> k_T_pm_raw;
    real<lower=0> theta_L_raw;
    real<lower=0,upper=1> theta_IPTG_raw;
    real<lower=0> n_IPTG_raw;
    real<lower=0> n_L_raw;

}

transformed parameters {
  // Introduction of the paramemters in an indexed object
  real theta[nParms];
  
  theta[1] = (k_in_IPTG_raw+2)*((0.4-0.004)/(4.00))+0.004;
  theta[2] = (k_out_IPTG_raw+2)*((0.4-0.004)/(4.00))+0.004;
  theta[3] = (k_in_aTc_raw+2)*((1-0.01)/(4.00))+0.01;
  theta[4] = (k_out_aTc_raw+2)*((1-0.01)/(4.00))+0.01;
  theta[5] = (k_L_pm0_raw+2)*((0.3-0.003)/(4.00))+0.003;
  theta[6] = (k_L_pm_raw+2)*((100.00-0.00)/(4.00))+1.00;
  theta[7] = (theta_T_raw+2)*((300.00-3.00)/(4.00))+3.00;
  theta[8] = (theta_aTc_raw+2)*((100.00-0.00)/(4.00))+0.00;
  theta[9] = (n_aTc_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[10] = (n_T_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[11] = (k_T_pm0_raw+2)*((1.00-0.01)/(4.00))+0.01;
  theta[12] = (k_T_pm_raw+2)*((10.00-0.1)/(4.00))+0.1;
  theta[13] = (theta_L_raw+2)*((300.00-3.00)/(4.00))+3.00;
  theta[14] = (theta_IPTG_raw+2)*((1.00-0.00)/(4.00))+0.00;
  theta[15] = (n_IPTG_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[16] = (n_L_raw+2)*((5.00-0.00)/(4.00))+0.00;

}

model {
  // Intermediate parameters
  real y_hat[tsl,Neq]; // Object that will include the solution for the ODEs for all the events
  real initialV[Neq]; // Initial value for the ODE each event
  int i; // Increasing index for the inputs
  vector[2] y_al; // Vector that will include the solutio of the algebraic solution for the steady state of the model
  real Y0[Neq]; // Initial values for the ODEs variables at the first event
  real ssv[tonil,Neq];
  
  // Priors definition 
  k_in_IPTG_raw ~ normal(0.202,0.099);
  k_out_IPTG_raw ~ normal(0.202,0.099);
  k_in_aTc_raw ~ normal(0.505,0.2475);
  k_out_aTc_raw ~ normal(0.505,0.2475);
  k_L_pm0_raw ~ normal(0.1515,0.07425);
  k_L_pm_raw ~ normal(50.5,24.75);
  theta_T_raw ~ normal(151.5,74.25);
  theta_aTc_raw ~ normal(50,25);
  n_aTc_raw ~ normal(2.5,1.25);
  n_T_raw ~ normal(2.5,1.25);
  k_T_pm0_raw ~ normal(0.505,0.2475);
  k_T_pm_raw ~ normal(5.05,2.475);
  theta_L_raw ~ normal(151.5,74.25);
  theta_IPTG_raw ~ normal(0.5,0.25);
  n_IPTG_raw ~ normal(2.5,1.25);
  n_L_raw ~ normal(2.5,1.25);

  // Calculation of initial guesses
  y_al = SteadyState(to_vector(ivss), to_vector(theta), pre, x_i); // Calculation of initial guesses for steady state
  Y0[1] = preIPTG;
  Y0[2] = preaTc;
  Y0[3] = y_al[1];
  Y0[4] = y_al[2];
  ssv = integrate_ode_bdf(Toogle_one, Y0,0,toni,theta,pre, x_i, 1e-9, 1e-9, 1e7); // ON incubation calculation for the steady state
  Y0 = ssv[tonil];

  // ODE simulation
  initialV = Y0;
  i = 1;
    // Loop (over the number of events) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    for (q in 1:Nevents){
      int itp = evnT[q];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q]+1):(evnT[q+1]+1)]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop
    
      // Calculation of the solution for the ODEs where for events that are not the firt one the time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){part1 = integrate_ode_bdf(Toogle_one, initialV,itp,ts[(evnT[q]+1):(evnT[q+1]+1)],theta,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);
      }
      else{part1 = integrate_ode_bdf(Toogle_one, initialV,(itp-1e-7),ts[(evnT[q]+1):(evnT[q+1]+1)],theta,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);
      }

      // Modification of the initial state values for the next event
      initialV = part1[lts];
      // Increase index for inputs
      i=i+2;
      // Introduction of the result of part1 into the object y_hat
      for (y in (itp+1):(itp+lts)){
        y_hat[(y),]=(part1)[(y-itp),];
      };
    };

  // Likelihood at each sampling time
  for (t in 1:stsl){
    RFPmean[t] ~ normal(y_hat[(sts[t]+1),3],RFPstd[t]);
    GFPmean[t] ~ normal(y_hat[(sts[t]+1),4],GFPstd[t]);
  }

}
