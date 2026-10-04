#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ReitsmaSubgroup(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(TP);
    DATA_VECTOR(FP);
    DATA_VECTOR(FN);
    DATA_VECTOR(TN);
    DATA_FACTOR(group);     // subgroup, 0, ..., G-1

    int G = 0;
    for(int i=0; i<group.size(); i++)
    {
        if(group(i) > G)
           G = group(i);
    }
    G++;

    int N = TP.size();
    
    if(FP.size() != N)
        error("FP and TP must have equal lengths");

    if(FN.size() != N)
        error("FN and TP must have equal lengths");

    if(TN.size() != N)
        error("TN and TP must have equal lengths");

    if(group.size() != N)
        error("group and TP must have equal lengths");

    /* ===== FIXED EFFECTS ===== */
    PARAMETER_VECTOR(mu_A);
    PARAMETER_VECTOR(mu_B);

    if(mu_A.size() != G)
      error("mu_A length must equal number of groups");

    if(mu_B.size() != G)
      error("mu_B length must equal number of groups");

    PARAMETER_VECTOR(log_sigma_A);
    PARAMETER_VECTOR(log_sigma_B);
    PARAMETER_VECTOR(theta_AB);

    if(log_sigma_A.size() != G)
      error("log_sigma_A length must equal number of groups");

    if(log_sigma_B.size() != G)
      error("log_sigma_B length must equal number of groups");

    if(theta_AB.size() != G)
      error("theta_AB length must equal number of groups");

    /* ===== RANDOM EFFECTS ===== */
    PARAMETER_VECTOR(sensu);
    PARAMETER_VECTOR(specu);

    if(sensu.size() != N)
        error("sensu length must equal number of studies");

    if(specu.size() != N)
        error("specu length must equal number of studies");

    
    /* ===== NEGATIVE LOGLIKELIHOOD ===== */
    Type nll = 0;

    /* ===== RANDOM EFFECTS VARIANCE-COVARIANCE MATRIX ===== */
    vector<matrix<Type>> Sigma(G);

    vector<Type> sigma_A(G);
    vector<Type> sigma_B(G);
    vector<Type> rho_AB(G);
    
    for(int g=0; g<G; g++)
    {
        sigma_A(g) = exp(log_sigma_A(g));
        sigma_B(g) = exp(log_sigma_B(g));
        rho_AB(g)  = Type(0.9999) * tanh(theta_AB(g));

        Sigma(g).resize(2,2);
  
        Sigma(g)(0,0) = sigma_A(g) * sigma_A(g);
        Sigma(g)(1,1) = sigma_B(g) * sigma_B(g);
        Sigma(g)(0,1) = rho_AB(g) * sigma_A(g) * sigma_B(g);
        Sigma(g)(1,0) = Sigma(g)(0,1);
    }

    for(int i=0;i<N;i++)
    {   
        int g = group(i);
    
        MVNORM_t<Type> neg_log_density(Sigma(g));
        
        vector<Type> ui(2);
        ui(0) = sensu(i);
        ui(1) = specu(i);
        
        nll += neg_log_density(ui);
    
       /* ===== LINEAR PREDICTORS ===== */
       Type eta_lsens = mu_A(g) + sensu(i);
       Type eta_lspec = mu_B(g) + specu(i);

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
    vector<Type> sigma2_A(G);
    vector<Type> sigma2_B(G);
    vector<Type> sigma_AB(G);
    
    for(int g=0; g<G; g++) 
    {
       sigma2_A(g) = sigma_A(g) * sigma_A(g);
       sigma2_B(g) = sigma_B(g) * sigma_B(g);
       sigma_AB(g) = rho_AB(g) * sigma_A(g) * sigma_B(g);
    }
    
    vector<Type> nu_A(G-1);
    vector<Type> nu_B(G-1); 
    if(G > 1)
    {  
       for(int g=1; g<G; g++)
       {
          nu_A(g-1) = mu_A(g) - mu_A(0);
          nu_B(g-1) = mu_B(g) - mu_B(0);
       }
    }

    REPORT(mu_A);
    REPORT(mu_B);
    REPORT(sigma2_A);
    REPORT(sigma2_B);
    REPORT(sigma_AB);
    REPORT(rho_AB);
    
    if(G > 1){
    REPORT(nu_A);
    REPORT(nu_B);  
    }

    ADREPORT(mu_A);
    ADREPORT(mu_B);
    ADREPORT(sigma2_A);
    ADREPORT(sigma2_B);
    ADREPORT(sigma_AB);
    ADREPORT(rho_AB);

    if(G > 1){
    ADREPORT(nu_A);
    ADREPORT(nu_B);  
    }

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this