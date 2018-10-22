# Function to calculate steady state for model 1 for initial guess of the initial values for each variable

ss1 <- function (v){
  
  # Parameters and input values have to be modified inside the function due to the non linear analytical solver used
  p<-c(2e-2, 2.75e-2, 1.11e-1, 3.2e-2, 8.3, 30, 11.65, 2, 2, 1.386e-1, 1.19e-1, 2.06, 39.94, 9.06e-2, 2, 2, 1.386e-1, 9.726e-1, 1.65e-1, 1.170, 1.65e-2)
  x_r<-c(0,0)
    u_IPTG = x_r[1];
    u_aTc = x_r[2];
  
    kml0 = p[1];
    kml = p[2];
    theta_T = p[3];
    theta_aTc = p[4];
    n_aTc = p[5];
    n_T = p[6];
    g_ml = p[7];
    kmt0 = p[8];
    kmt = p[9];
    theta_L = p[10];
    theta_IPTG = p[11];
    n_IPTG = p[12];
    n_L = p[13];
    g_mt = p[14];
    k_pl = p[15];
    g_pl = p[16];
    k_pt = p[17];
    g_pt = p[18];
    
    eq <- numeric(4)

    eq[1] = p[1]+p[2]/(1+(v[4]/p[3]*1/(1+(x_r[2]/p[4])^p[5]))^p[6])-p[7]*v[1];
    
    eq[2] = p[8]+p[9]/(1+(v[3]/p[10]*1/(1+(x_r[1]/p[11])^p[12]))^p[13])-p[14]*v[2];
    eq[3] = p[15]*v[1]-p[16]*v[3];
    eq[4] = p[17]*v[2]-p[18]*v[4];
    return(eq);
 
}

# nleqslv need to be installed, not installed with R
nleqslv(c(0,0,0,0),ss1)->q
transform(q)[,1]->w


