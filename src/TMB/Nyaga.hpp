#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type Nyaga(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(TP);
    DATA_VECTOR(FP);
    DATA_VECTOR(FN);
    DATA_VECTOR(TN);
    DATA_FACTOR(study);  // 0, ..., I-1
    DATA_FACTOR(test);   // 0, ..., K-1

    int N = TP.size();

    if(FP.size() != N)
        error("FP and TP must have equal lengths");

    if(FN.size() != N)
        error("FN and TP must have equal lengths");

    if(TN.size() != N)
        error("TN and TP must have equal lengths");

    if(study.size() != N)
        error("study and TP must have equal lengths");

    if(test.size() != N)
        error("test and TP must have equal lengths");

    int I = 0;
    for(int r = 0; r < N; r++)
        if(study(r) > I)
            I = study(r);
    I++;

    int K = 0;
    for(int r = 0; r < N; r++)
        if(test(r) > K)
            K = test(r);
    K++;

    /* ===== FIXED EFFECTS ===== */
    // Test-specific mean logit sensitivity and specificity
    PARAMETER_VECTOR(mu_A);
    PARAMETER_VECTOR(mu_B);

    if(mu_A.size() != K)
        error("mu_A length must equal number of tests");

    if(mu_B.size() != K)
        error("mu_B length must equal number of tests");

    /* ===== STUDY-LEVEL HETEROGENEITY ===== */
    PARAMETER(log_sigma_A);
    PARAMETER(log_sigma_B);
    PARAMETER(theta_AB);

    Type sigma_A = exp(log_sigma_A);
    Type sigma_B = exp(log_sigma_B);
    Type rho_AB  = Type(0.9999) * tanh(theta_AB);

    matrix<Type> Sigma(2, 2);

    Sigma(0, 0) = sigma_A * sigma_A;
    Sigma(1, 1) = sigma_B * sigma_B;
    Sigma(0, 1) = rho_AB * sigma_A * sigma_B;
    Sigma(1, 0) = Sigma(0, 1);

    MVNORM_t<Type> study_density(Sigma);

    // One pair per study, not one pair per arm
    PARAMETER_VECTOR(sensu);
    PARAMETER_VECTOR(specu);

    if(sensu.size() != I)
        error("sensu length must equal number of studies");

    if(specu.size() != I)
        error("specu length must equal number of studies");

    /* ===== ARM-LEVEL RESIDUAL HETEROGENEITY ===== */
    PARAMETER(log_tau_A);
    PARAMETER(log_tau_B);

    Type tau_A = exp(log_tau_A);
    Type tau_B = exp(log_tau_B);

    // One residual deviation per observed study-test arm
    PARAMETER_VECTOR(sens_delta);
    PARAMETER_VECTOR(spec_delta);

    if(sens_delta.size() != N)
        error("sens_delta length must equal number of arms");

    if(spec_delta.size() != N)
        error("spec_delta length must equal number of arms");

    /* ===== NEGATIVE LOG-LIKELIHOOD ===== */
    Type nll = Type(0.0);

    // Study-level correlated random effects
    for(int i = 0; i < I; i++)
    {
        vector<Type> ui(2);

        ui(0) = sensu(i);
        ui(1) = specu(i);

        nll += study_density(ui);
    }

    // Arm-level independent residual random effects
    for(int r = 0; r < N; r++)
    {
        nll -= dnorm(sens_delta(r), Type(0.0), tau_A, true);
        nll -= dnorm(spec_delta(r), Type(0.0), tau_B, true);
    }

    // Conditional binomial likelihood
    for(int r = 0; r < N; r++)
    {
        int i = study(r);
        int k = test(r);

        Type eta_lsens = mu_A(k) + sensu(i) + sens_delta(r);
        Type eta_lspec = mu_B(k) + specu(i) + spec_delta(r);

        Type sens = invlogit(eta_lsens);
        Type spec = invlogit(eta_lspec);

        Type n_sens = TP(r) + FN(r);
        Type n_spec = TN(r) + FP(r);

        nll -= dbinom(TP(r), n_sens, sens, true);
        nll -= dbinom(TN(r), n_spec, spec, true);
    }

    /* ===== DERIVED QUANTITIES ===== */
    Type sigma2_A  = sigma_A * sigma_A;
    Type sigma2_B  = sigma_B * sigma_B;
    Type sigma_AB  = rho_AB * sigma_A * sigma_B;

    Type tau2_A = tau_A * tau_A;
    Type tau2_B = tau_B * tau_B;

    // Total arm-level marginal variance
    Type total_var_A = sigma2_A + tau2_A;
    Type total_var_B = sigma2_B + tau2_B;

    // Correlation between two different sensitivity or specificity arms observed in the same study
    Type icc_A = sigma2_A / total_var_A;
    Type icc_B = sigma2_B / total_var_B;

    // Cross-outcome correlation for two arm-level latent predictors
    Type total_rho_AB = sigma_AB / sqrt(total_var_A * total_var_B);

    vector<Type> nu_A(K-1);
    vector<Type> nu_B(K-1);
    if(K > 1)
    {  
       for(int k=1; k<K; k++)
       {
          nu_A(k-1) = mu_A(k) - mu_A(0);
          nu_B(k-1) = mu_B(k) - mu_B(0);
       }
    }

    /* ===== REPORTS ===== */

    REPORT(mu_A);
    REPORT(mu_B);

    REPORT(sigma2_A);
    REPORT(sigma2_B);
    REPORT(sigma_AB);
    REPORT(rho_AB);

    REPORT(tau2_A);
    REPORT(tau2_B);

    REPORT(total_var_A);
    REPORT(total_var_B);
    REPORT(icc_A);
    REPORT(icc_B);
    REPORT(total_rho_AB);
    
    if(K > 1){
    REPORT(nu_A);
    REPORT(nu_B);
    }

    ADREPORT(mu_A);
    ADREPORT(mu_B);

    ADREPORT(sigma2_A);
    ADREPORT(sigma2_B);
    ADREPORT(sigma_AB);
    ADREPORT(rho_AB);

    ADREPORT(tau2_A);
    ADREPORT(tau2_B);

    ADREPORT(total_var_A);
    ADREPORT(total_var_B);
    ADREPORT(icc_A);
    ADREPORT(icc_B);
    ADREPORT(total_rho_AB);
    
    if(K > 1){
    ADREPORT(nu_A);
    ADREPORT(nu_B);
    }

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this