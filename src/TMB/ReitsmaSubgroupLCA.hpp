#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ReitsmaSubgroupLCA(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y11);
    DATA_VECTOR(y10);
    DATA_VECTOR(y01);
    DATA_VECTOR(y00);
    DATA_FACTOR(group);     // subgroup, 0, ..., G-1

    int G = 0;
    for(int i=0; i<group.size(); i++)
    {
        if(group(i) > G)
           G = group(i);
    }
    G++;

    int N = y11.size();
    
    if(y10.size() != N)
        error("y10 and y11 must have equal lengths");

    if(y01.size() != N)
        error("y01 and y11 must have equal lengths");

    if(y00.size() != N)
        error("y00 and y11 must have equal lengths");

    if(group.size() != N)
        error("group and y11 must have equal lengths");

    /* ===== FIXED EFFECTS ===== */
    PARAMETER_VECTOR(mu_prev);
    PARAMETER_VECTOR(mu_A_index);
    PARAMETER_VECTOR(mu_B_index);

    if(mu_prev.size() != G)
      error("mu_prev length must equal number of groups");

    if(mu_A_index.size() != G)
      error("mu_A_index length must equal number of groups");

    if(mu_B_index.size() != G)
      error("mu_B_index length must equal number of groups");

    PARAMETER(mu_A_ref);
    PARAMETER(mu_B_ref);

    PARAMETER_VECTOR(log_sigma_prev);
    PARAMETER_VECTOR(log_sigma_A_index);
    PARAMETER_VECTOR(log_sigma_B_index);
    PARAMETER_VECTOR(theta_AB_index);

    if(log_sigma_prev.size() != G)
      error("log_sigma_prev length must equal number of groups");

    if(log_sigma_A_index.size() != G)
      error("log_sigma_A_index length must equal number of groups");

    if(log_sigma_B_index.size() != G)
      error("log_sigma_B_index length must equal number of groups");

    if(theta_AB_index.size() != G)
      error("theta_AB_index length must equal number of groups");

    /* ===== RANDOM EFFECTS ===== */
    PARAMETER_VECTOR(prevu);
    PARAMETER_VECTOR(sensu);
    PARAMETER_VECTOR(specu);

    if(prevu.size() != N)
        error("prevu length must equal number of studies");

    if(sensu.size() != N)
        error("sensu length must equal number of studies");

    if(specu.size() != N)
        error("specu length must equal number of studies");

    
    /* ===== NEGATIVE LOGLIKELIHOOD ===== */
    Type nll = 0;
    Type eps = Type(1e-12);

    /* ===== LINEAR PREDICTORS AND TRANSFORM ===== */
    Type eta_lsensref   = mu_A_ref;
    Type eta_lspecref   = mu_B_ref;
    Type sensref   = invlogit(eta_lsensref);
    Type specref   = invlogit(eta_lspecref);

    /* ===== RANDOM EFFECTS VARIANCE-COVARIANCE MATRIX ===== */
    vector< matrix<Type> > Sigma(G);

    for(int g=0; g<G; g++)
    {
        Type sigma_prev    = exp(log_sigma_prev(g));
        Type sigma_A_index = exp(log_sigma_A_index(g));
        Type sigma_B_index = exp(log_sigma_B_index(g));
        Type rho_AB        = Type(0.9999) * tanh(theta_AB_index(g));

        Sigma(g).resize(3,3);
        Sigma(g).setZero();

        Sigma(g)(0,0) = sigma_prev * sigma_prev;
        Sigma(g)(1,1) = sigma_A_index * sigma_A_index;
        Sigma(g)(2,2) = sigma_B_index * sigma_B_index;
        Sigma(g)(1,2) = rho_AB * sigma_A_index * sigma_B_index;
        Sigma(g)(2,1) = Sigma(g)(1,2);
    }

    for(int i=0;i<N;i++)
    {   
        int g = group(i);
    
        MVNORM_t<Type> neg_log_density(Sigma(g));
        
        vector<Type> ui(3);
        ui(0) = prevu(i);
        ui(1) = sensu(i);
        ui(2) = specu(i);
        
        nll += neg_log_density(ui);
    
       /* ===== LINEAR PREDICTORS ===== */
       Type eta_prev       = mu_prev(g)    + prevu(i);
       Type eta_lsensindex = mu_A_index(g) + sensu(i);
       Type eta_lspecindex = mu_B_index(g) + specu(i);

       /* ===== TRANSFORM ===== */
       Type prev      = invlogit(eta_prev);
       Type sensindex = invlogit(eta_lsensindex);
       Type specindex = invlogit(eta_lspecindex);

       /* ===== LATENT CLASS PROBABILITIES ===== */
       Type p11 = prev * (sensindex * sensref) + (Type(1.0) - prev) * (Type(1.0) - specindex) * (Type(1.0) - specref);
       Type p10 = prev * (sensindex * (Type(1.0) - sensref)) + (Type(1.0) - prev) * (Type(1.0) - specindex) * specref;
       Type p01 = prev * (Type(1.0) - sensindex) * sensref   + (Type(1.0) - prev) * specindex * (Type(1.0) - specref);
       Type p00 = prev * (Type(1.0) - sensindex) * (Type(1.0) - sensref) + (Type(1.0) - prev) * (specindex * specref);

       p11 = CppAD::CondExpLt(p11, eps, eps, p11);
       p10 = CppAD::CondExpLt(p10, eps, eps, p10);
       p01 = CppAD::CondExpLt(p01, eps, eps, p01);
       p00 = CppAD::CondExpLt(p00, eps, eps, p00);

       /* ===== LIKELIHOOD CONTRIBUTIONS ===== */
       nll -= y11(i) * log(p11) + y10(i) * log(p10) + y01(i) * log(p01) + y00(i) * log(p00);

    }

    /* ===== REPORTS ===== */

    vector<Type> sigma2_prev(G);
    vector<Type> sigma2_A_index(G);
    vector<Type> sigma2_B_index(G);
    vector<Type> sigma_AB_index(G);
    vector<Type> rho_AB_index(G);

    for(int g=0; g<G; g++) 
    {
       sigma2_prev(g)    = exp(Type(2.0)*log_sigma_prev(g));
       sigma2_A_index(g) = exp(Type(2.0)*log_sigma_A_index(g));
       sigma2_B_index(g) = exp(Type(2.0)*log_sigma_B_index(g));
       rho_AB_index(g)   = Type(0.9999) * tanh(theta_AB_index(g));
       sigma_AB_index(g) = rho_AB_index(g) * exp(log_sigma_A_index(g)) * exp(log_sigma_B_index(g));
    }
    
    vector<Type> nu_prev(G-1);
    vector<Type> nu_A_index(G-1);
    vector<Type> nu_B_index(G-1); 
    if(G > 1)
    {  
       for(int g=1; g<G; g++)
       {
          nu_prev(g-1)    = mu_prev(g) - mu_prev(0);
          nu_A_index(g-1) = mu_A_index(g) - mu_A_index(0);
          nu_B_index(g-1) = mu_B_index(g) - mu_B_index(0);
       }
    }

    REPORT(mu_prev);
    REPORT(mu_A_index);
    REPORT(mu_B_index);

    REPORT(sigma2_prev);
    REPORT(sigma2_A_index);
    REPORT(sigma2_B_index);
    REPORT(sigma_AB_index);
    REPORT(rho_AB_index);
    
    REPORT(mu_A_ref);
    REPORT(mu_B_ref);

    if(G > 1){
    REPORT(nu_prev);
    REPORT(nu_A_index);
    REPORT(nu_B_index);  
    }

    ADREPORT(mu_prev);
    ADREPORT(mu_A_index);
    ADREPORT(mu_B_index);

    ADREPORT(sigma2_prev);
    ADREPORT(sigma2_A_index);
    ADREPORT(sigma2_B_index);
    ADREPORT(sigma_AB_index);
    ADREPORT(rho_AB_index);
    
    ADREPORT(mu_A_ref);
    ADREPORT(mu_B_ref);

    if(G > 1){
    ADREPORT(nu_prev);
    ADREPORT(nu_A_index);
    ADREPORT(nu_B_index);  
    }

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this