#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type Reitsma(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(TP);
    DATA_VECTOR(FP);
    DATA_VECTOR(FN);
    DATA_VECTOR(TN);

    int N = TP.size();
    
    if(FP.size() != N)
        error("FP and TP must have equal lengths");

    if(FN.size() != N)
        error("FN and TP must have equal lengths");

    if(TN.size() != N)
        error("TN and TP must have equal lengths");

    /* ===== FIXED EFFECTS ===== */
    PARAMETER(mu_A);
    PARAMETER(mu_B);

    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    PARAMETER(theta_AB);
    
    Type sigma_A = exp(log_sigma_A);
    Type sigma_B = exp(log_sigma_B);
    Type rho_AB  = Type(0.9999) * tanh(theta_AB);

    /* ===== RANDOM EFFECTS ===== */
    PARAMETER_VECTOR(sensu);
    PARAMETER_VECTOR(specu);

    if(sensu.size() != N)
        error("sensu length must equal number of studies");

    if(specu.size() != N)
        error("specu length must equal number of studies");

    matrix<Type> Sigma(2,2);
    Sigma(0,0) = sigma_A * sigma_A;
    Sigma(1,1) = sigma_B * sigma_B;
    Sigma(0,1) = rho_AB * sigma_A * sigma_B;
    Sigma(1,0) = Sigma(0,1);
    
    MVNORM_t<Type> neg_log_density(Sigma);
    
    /* ===== NEGATIVE LOGLIKELIHOOD ===== */
    Type nll = Type(0.0);

    for(int i=0;i<N;i++)
    {
        vector<Type> ui(2);

        ui(0) = sensu(i);
        ui(1) = specu(i);
        
        nll += neg_log_density(ui);
    
       /* ===== LINEAR PREDICTORS ===== */
       Type eta_lsens = mu_A + sensu(i);
       Type eta_lspec = mu_B + specu(i);

       /* ===== TRANSFORM ===== */
       Type sens = invlogit(eta_lsens);
       Type spec = invlogit(eta_lspec);
   
       /* ===== LIKELIHOOD CONTRIBUTIONS ===== */
       Type n_sens = TP(i) + FN(i);
       Type n_spec = TN(i) + FP(i);
       nll -= dbinom(TP(i),n_sens,sens,true);
       nll -= dbinom(TN(i),n_spec,spec,true);

    }

    /* ===== REPORTS ===== */
    Type sigma2_A = sigma_A * sigma_A;
    Type sigma2_B = sigma_B * sigma_B;
    Type sigma_AB = rho_AB * sigma_A * sigma_B;

    REPORT(mu_A);
    REPORT(mu_B);
    REPORT(sigma2_A);
    REPORT(sigma2_B);
    REPORT(sigma_AB);
    REPORT(rho_AB);
    
    ADREPORT(mu_A);
    ADREPORT(mu_B);
    ADREPORT(sigma2_A);
    ADREPORT(sigma2_B);
    ADREPORT(sigma_AB);
    ADREPORT(rho_AB);

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this