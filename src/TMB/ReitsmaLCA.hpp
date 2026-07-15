#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ReitsmaLCA(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y11);
    DATA_VECTOR(y10);
    DATA_VECTOR(y01);
    DATA_VECTOR(y00);

    int N = y11.size();
    
    if(y10.size() != N)
        error("y10 and y11 must have equal lengths");

    if(y01.size() != N)
        error("y01 and y11 must have equal lengths");

    if(y00.size() != N)
        error("y00 and y11 must have equal lengths");

    /* ===== FIXED EFFECTS ===== */
    PARAMETER(mu_prev);
    PARAMETER(mu_A_index);
    PARAMETER(mu_B_index);
    PARAMETER(mu_A_ref);
    PARAMETER(mu_B_ref);

    PARAMETER(log_sigma_prev);
    PARAMETER(log_sigma_A_index);
    PARAMETER(log_sigma_B_index);
    PARAMETER(theta_AB_index);
    
    Type sigma_prev    = exp(log_sigma_prev);
    Type sigma_A_index = exp(log_sigma_A_index);
    Type sigma_B_index = exp(log_sigma_B_index);
    Type rho_AB_index  = Type(0.9999) * tanh(theta_AB_index);

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

    matrix<Type> Sigma(3,3);
    Sigma.setZero();
    Sigma(0,0) = sigma_prev * sigma_prev;
    Sigma(1,1) = sigma_A_index * sigma_A_index;
    Sigma(2,2) = sigma_B_index * sigma_B_index;
    Sigma(1,2) = rho_AB_index * sigma_A_index * sigma_B_index;
    Sigma(2,1) = Sigma(1,2);
    
    MVNORM_t<Type> neg_log_density(Sigma);
    
    /* ===== NEGATIVE LOGLIKELIHOOD ===== */
    Type nll = 0;
    Type eps = Type(1e-12);

    /* ===== LINEAR PREDICTORS AND TRANSFORM ===== */
    Type eta_lsensref   = mu_A_ref;
    Type eta_lspecref   = mu_B_ref;
    Type sensref   = invlogit(eta_lsensref);
    Type specref   = invlogit(eta_lspecref);

    for(int i=0;i<N;i++)
    {
        vector<Type> ui(3);

        ui(0) = prevu(i);
        ui(1) = sensu(i);
        ui(2) = specu(i);
        
        nll += neg_log_density(ui);
    
       /* ===== LINEAR PREDICTORS ===== */
       Type eta_prev       = mu_prev    + prevu(i);
       Type eta_lsensindex = mu_A_index + sensu(i);
       Type eta_lspecindex = mu_B_index + specu(i);

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
    Type sigma2_prev    = sigma_prev * sigma_prev;
    Type sigma2_A_index = sigma_A_index * sigma_A_index;
    Type sigma2_B_index = sigma_B_index * sigma_B_index;
    Type sigma_AB_index = rho_AB_index * sigma_A_index * sigma_B_index;

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

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this