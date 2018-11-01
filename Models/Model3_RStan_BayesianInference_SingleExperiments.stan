// RStan model containing the inducer exchange model for the Toggle Switch proposed by cSynBioSys

functions{
  // Function containing the ODE to be used for the inference
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
     
    // INPUTS
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    //PARAMETERS
    real k_IPTG = p[1];
    real k_aTc = p[2];
    real k_L_pm0 = p[3];
    real k_L_pm = p[4];
    real theta_T = p[5];
    real theta_aTc = p[6];
    real n_aTc = p[7];
    real n_T = p[8];
    real k_T_pm0 = p[9];
    real k_T_pm = p[10];
    real theta_L = p[11];
    real theta_IPTG = p[12];
    real n_IPTG = p[13];
    real n_L = p[14];
    
    // EQUATIONS
    real dInd_dt[4];

    dInd_dt[1] = k_IPTG*(x_r[1]-y[1])-0.0165*y[1];
    dInd_dt[2] = k_aTc*(x_r[2]-y[2])-0.0165*y[2];

    dInd_dt[3] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+(y[4]/theta_T*1/(1+(y[2]/theta_aTc)^n_aTc))^n_T))))-0.0165*y[3];
    dInd_dt[4] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+(y[3]/theta_L*1/(1+(y[1]/theta_IPTG)^n_IPTG))^n_L))))-0.0165*y[4];

    // RESULTS
    return dInd_dt;

  }
  
  // Function type vector containing the equations where the root needs to be calculated for the steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    vector[4] alpha;
    // Parameters

    real k_IPTG = p[1];
    real k_aTc = p[2];
    real k_L_pm0 = p[3];
    real k_L_pm = p[4];
    real theta_T = p[5];
    real theta_aTc = p[6];
    real n_aTc = p[7];
    real n_T = p[8];
    real k_T_pm0 = p[9];
    real k_T_pm = p[10];
    real theta_L = p[11];
    real theta_IPTG = p[12];
    real n_IPTG = p[13];
    real n_L = p[14];
    
    // Equations
    alpha[1] = (k_IPTG*init[1])/(k_IPTG+0.0165);
    alpha[2] = (k_aTc*init[2])/(k_aTc+0.0165);
    
    alpha[3] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+((init[4]/theta_T)*1/(1+(alpha[2]/theta_aTc)^n_aTc))^n_T))))/0.0165;
    alpha[4] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+((init[3]/theta_L)*1/(1+(alpha[1]/theta_IPTG)^n_IPTG))^n_L))))/0.0165;
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
  int nParms = 14; // Number of parameters of the model
  int Neq = 4; // Total number of equations of the model
  int Nevents = Nsp-1; // Number of events
  int x_i[0]; // Empty x_i object (needs to be deffined)
  real x_r[(Nevents)*2]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq]; // Initial experimental values for the calculation of the steady state ordered as IPTG, aTc
  real pre[2]; // Input values during the ON incubation ordered as IPTG, aTc
  
  ivss[1] = preIPTG;
  ivss[2] = preaTc;
  ivss[3] = RFPmean[1];
  ivss[4] = GFPmean[1];
  
  pre[1] = preIPTG;
  pre[2] = preaTc;
  
}

parameters {
    // Parameters to be infered in the model
    real<lower=0> k_IPTG_raw;
    real<lower=0> k_aTc_raw;
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
  
  theta[1] = (k_IPTG_raw+2)*((0.4-0.004)/(4.00))+0.004;
  theta[2] = (k_aTc_raw+2)*((1-0.01)/(4.00))+0.01;
  theta[3] = (k_L_pm0_raw+2)*((0.3-0.003)/(4.00))+0.003;
  theta[4] = (k_L_pm_raw+2)*((100.00-0.00)/(4.00))+1.00;
  theta[5] = (theta_T_raw+2)*((300.00-3.00)/(4.00))+3.00;
  theta[6] = (theta_aTc_raw+2)*((100.00-0.00)/(4.00))+0.00;
  theta[7] = (n_aTc_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[8] = (n_T_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[9] = (k_T_pm0_raw+2)*((1.00-0.01)/(4.00))+0.01;
  theta[10] = (k_T_pm_raw+2)*((10.00-0.1)/(4.00))+0.1;
  theta[11] = (theta_L_raw+2)*((300.00-3.00)/(4.00))+3.00;
  theta[12] = (theta_IPTG_raw+2)*((1.00-0.00)/(4.00))+0.00;
  theta[13] = (n_IPTG_raw+2)*((5.00-0.00)/(4.00))+0.00;
  theta[14] = (n_L_raw+2)*((5.00-0.00)/(4.00))+0.00;

}

model {
  // Intermediate parameters
  real y_hat[tsl,Neq]; // Object that will include the solution for the ODEs for all the events
  real initialV[Neq]; // Initial value for the ODE each event
  int i; // Increasing index for the inputs
  vector[4] y_al; // Vector that will include the solutio of the algebraic solution for the steady state of the model
  real Y0[Neq]; // Initial values for the ODEs variables at the first event
  real ssv[tonil,Neq];
  
  // Priors definition 
  k_IPTG_raw ~ normal(0.202,0.099);
  k_aTc_raw ~ normal(0.505,0.2475);
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
  Y0[1] = y_al[1];
  Y0[2] = y_al[2];
  Y0[3] = y_al[3];
  Y0[4] = y_al[4];
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
      if (q == 1){part1 = integrate_ode_bdf(Toogle_one, initialV,itp,ts[(evnT[q]+1):(evnT[q+1]+1)],theta,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);}
      else{part1 = integrate_ode_bdf(Toogle_one, initialV,(itp-1e-7),ts[(evnT[q]+1):(evnT[q+1]+1)],theta,to_array_1d(inputs[i:(i+1)]), x_i, 1e-9, 1e-9, 1e7);}

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
