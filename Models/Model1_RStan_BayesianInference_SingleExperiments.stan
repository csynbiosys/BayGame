// RStan model containing the  first inducer exchange model for the Toggle Switch contained in the Lugagne paper

functions{
  // Function containing the ODE to be used for the inference
  real[] Toogle_one(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    // Inputs
    real u_IPTG = x_r[1];
    real u_aTc = x_r[2];
    // Parameters
    real k_IPTG = p[1];
    real k_L_pm0 = p[2];
    real k_L_pm = p[3];
    real theta_T = p[4];
    real theta_aTc = p[5];
    real n_aTc = p[6];
    real n_T = p[7];
    real k_T_pm0 = p[8];
    real k_T_pm = p[9];
    real theta_L = p[10];
    real theta_IPTG = p[11];
    real n_IPTG = p[12];
    real n_L = p[13];
    
    //Equations
    real dInd_dt[3];
    dInd_dt[1] = k_IPTG*(x_r[1]-y[1]);
    dInd_dt[2] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+(y[3]/theta_T*1/(1+(x_r[2]/theta_aTc)^n_aTc))^n_T))))-0.0165*y[2];
    dInd_dt[3] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+(y[2]/theta_L*1/(1+(y[1]/theta_IPTG)^n_IPTG))^n_L))))-0.0165*y[3];
    //RESULTS
    return dInd_dt;
  }
  
  // Function type vector containing the equations where the root needs to be calculated for the steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    // Parameters
    real k_IPTG = p[1];
    real k_L_pm0 = p[2];
    real k_L_pm = p[3];
    real theta_T = p[4];
    real theta_aTc = p[5];
    real n_aTc = p[6];
    real n_T = p[7];
    real k_T_pm0 = p[8];
    real k_T_pm = p[9];
    real theta_L = p[10];
    real theta_IPTG = p[11];
    real n_IPTG = p[12];
    real n_L = p[13];
    // Equations
    vector[2] alpha;
    alpha[1] = ((1/0.1386)*(k_L_pm0+(k_L_pm/(1+((init[2]/theta_T)*1/(1+(x_r[2]/theta_aTc)^n_aTc))^n_T))))-0.0165*init[1];
    alpha[2] = ((1/0.1386)*(k_T_pm0+(k_T_pm/(1+((init[1]/theta_L)*1/(1+(x_r[1]/theta_IPTG)^n_IPTG))^n_L))))-0.0165*init[2];
    // Results
    return alpha;
  }
  
}

data {
  // INPUTS
  int<lower = 0> tsl; // length of measurements/time series
  real ts[tsl]; // time series
  int tsmax; // max(ts)
  int preIPTG; // level of inputs during the ON
  int preaTc;
  int Nsp; // Number of events
  int evnT[Nsp]; // Event starting times
  real inputs[(Nsp-1)*2]; // Number of total inputs
  real IPTG[Nsp-1]; // levels of inputs at each timepoints
  real aTc[Nsp-1];
  
  // OBSERVABLES
  int<lower=0> stsl; // Length of the sampling time series
  int sts[stsl]; // sampling times
  real GFPmean[stsl]; // estimated observables for tetR+GFP at each sampling time
  real RFPmean[stsl]; // estimated observables for LacI+RFP at each sampling time
  real<lower=0> GFPstd[stsl]; // standard error for tetR+GFP at each sampling time
  real<lower=0> RFPstd[stsl]; // standard error for LacI+RFP at each sampling time

}

transformed data {
  int nParms = 13; // Number of parameters of the model
  int Neq = 3; // Total number of equations of the model
  int Nevents = num_elements(evnT)-1; // Total number of events
  int x_i[0]; // Empty x_i object (needs to be deffined)
  real x_r[(Nsp-1)*2]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq-1]; // Initial experimental values for the calculation of the steady state ordered as IPTG, aTc
  real pre[2]; // Input values during the ON incubation ordered as IPTG, aTc
  
  ivss[1] = RFPmean[1];
  ivss[2] = GFPmean[1];
  
  pre[1] = preIPTG;
  pre[2] = preaTc;
}

parameters {
    // Parameters to be infered in the model
    real<lower=0> k_IPTG;
    real<lower=0> k_L_pm0;
    real<lower=0> k_L_pm;
    real<lower=0> theta_T;
    real<lower=0> theta_aTc;
    real<lower=0> n_aTc;
    real<lower=0> n_T;
    real<lower=0> k_T_pm0;
    real<lower=0> k_T_pm;
    real<lower=0> theta_L;
    real<lower=0> theta_IPTG;
    real<lower=0> n_IPTG;
    real<lower=0> n_L;
    
}

transformed parameters {
  // Introduction of the paramemters in an indexed object
  real theta[nParms];
  theta[1] = k_IPTG;
  theta[2] = k_L_pm0;
  theta[3] = k_L_pm;
  theta[4] = theta_T;
  theta[5] = theta_aTc;
  theta[6] = n_aTc;
  theta[7] = n_T;
  theta[8] = k_T_pm0;
  theta[9] = k_T_pm;
  theta[10] = theta_L;
  theta[11] = theta_IPTG;
  theta[12] = n_IPTG;
  theta[13] = n_L;
  
}

model {
  // Intermediate parameters
  real y_hat[tsmax,Neq]; // Object that will include the solution for the ODEs for all the events
  real initialV[Neq]; // Initial value for the ODE each event
  int i; // Increasing index for the inputs
  vector[2] y_al; // Vector that will include the solutio of the algebraic solution for the steady state of the model
  real Y0[Neq]; // Initial values for the ODEs variables at the first event
  
  // Priors definition (test)
  k_IPTG ~ normal(0.04,0.00666666666666667);
  k_L_pm0 ~ normal(0.2011527,0.03352545);
  k_L_pm ~ normal(8.594406,1.432401);
  theta_T ~ normal(76.4,12.7333333333333);
  theta_aTc ~ normal(35.98,5.99666666666667);
  n_aTc ~ normal(2,0.333333333333333);
  n_T ~ normal(2.152,0.358666666666667);
  k_T_pm0 ~ normal(0.16889674,0.0281494566666667);
  k_T_pm ~ normal(2.577039,0.4295065);
  theta_L ~ normal(124.9,20.8166666666667);
  theta_IPTG ~ normal(0.2926,0.0487666666666667);
  n_IPTG ~ normal(2,0.333333333333333);
  n_L ~ normal(2,0.333333333333333);
  
  // Calculation of initial guesses
  y_al = algebra_solver(SteadyState, to_vector(ivss), to_vector(theta), pre, x_i);
  Y0[1] = preIPTG;
  Y0[2] = y_al[1];
  Y0[3] = y_al[2];
  
  // ODE simulation
  initialV = Y0;
  i = 1;
    // Loop (over the number of events) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    for (q in 1:Nevents){
      int itp = evnT[q];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q]+1):evnT[q+1]]);  // Length of the time series for each event
      real Tevent[lts] = ts[(evnT[q]+1):evnT[q+1]];  // Time series of each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop
      
      // Calculation of the solution for the ODEs where for events that are not the firt one the time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){part1 = integrate_ode_rk45(Toogle_one, initialV,itp,ts[(evnT[q]+1):evnT[q+1]],theta,to_array_1d(inputs[i:(i+1)]), x_i);}
      else{part1 = integrate_ode_rk45(Toogle_one, initialV,(itp-1),ts[(evnT[q]+1):evnT[q+1]],theta,to_array_1d(inputs[i:(i+1)]), x_i);}

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
    RFPmean[t] ~ normal(y_hat[(sts[t]+1),2],RFPstd[t]);
    GFPmean[t] ~ normal(y_hat[(sts[t]+1),3],GFPstd[t]);
  }
}
