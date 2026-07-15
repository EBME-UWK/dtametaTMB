#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type RutterGatsonisLCA(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y11);
    DATA_VECTOR(y10);
    DATA_VECTOR(y01);
    DATA_VECTOR(y00);
    DATA_VECTOR(spec);

    int N     = y11.size();
    int nspec = spec.size();
    
    if(y10.size() != N)
        error("y10 and y11 must have equal lengths");

    if(y01.size() != N)
        error("y01 and y11 must have equal lengths");

    if(y00.size() != N)
        error("y00 and y11 must have equal lengths");
        
    /* =====  VALIDATE SPEC GRID ===== */
    for(int i=0; i<nspec; i++)
    {
        double sp = asDouble(spec(i));
        if(!R_finite(sp)) error("spec contains non-finite values");
        if(sp <= 0.0 || sp >= 1.0) error("All values in spec must satisfy 0 < spec < 1");
    }
    
    /* ===== FIXED EFFECTS ===== */
    PARAMETER(mu_prev); // logit prevalence
    PARAMETER(Lambda);
    PARAMETER(Theta);
    PARAMETER(beta);
    PARAMETER(log_sigma_prev);
    PARAMETER(log_sigma_alpha);
    PARAMETER(log_sigma_theta);

    Type sigma_alpha  = exp(log_sigma_alpha);
    Type sigma_theta  = exp(log_sigma_theta);
    Type sigma_prev   = exp(log_sigma_prev); 
    
    /* ===== REFERENCE TEST ===== */
    PARAMETER(mu_A_ref);
    PARAMETER(mu_B_ref);
    
    Type sensref = invlogit(mu_A_ref);
    Type specref = invlogit(mu_B_ref);
    
    /* =====  RANDOM EFFECTS ===== */
    PARAMETER_VECTOR(prevu);
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(theta);
    
    if(prevu.size() != N)
        error("prevu length must equal number of studies");

    if(theta.size() != N)
        error("theta length must equal number of studies");

    if(alpha.size() != N)
        error("alpha length must equal number of studies");
    
    /* =====  NEGATIVE LOG-LIKELIHOOD ===== */
    Type nll = 0.0;
    Type eps = Type(1e-12);

    /* =====  RANDOM EFFECTS (VECTORIZED) ===== */
    nll -= sum(dnorm(alpha, Type(0.0), sigma_alpha, true));
    nll -= sum(dnorm(theta, Type(0.0), sigma_theta, true));
    nll -= sum(dnorm(prevu, Type(0.0), sigma_prev,  true));

    for(int i = 0; i < N; i++)
    {   
        /* ===== LINEAR PREDICTORS AND TRANSFORM ===== */
        Type prev = invlogit( mu_prev + prevu(i) );
        Type etaD = ( Theta + theta(i) + (Lambda + alpha(i))/Type(2.0) ) * exp(-beta/Type(2.0));
        Type etaH = ( Theta + theta(i) - (Lambda + alpha(i))/Type(2.0) ) * exp( beta/Type(2.0));
        Type sensindex = invlogit( etaD);
        Type specindex = invlogit(-etaH);
        
        /* ===== LATENT CLASS PROBABILITIES ===== */
        Type p11 = prev * (sensindex * sensref) + (Type(1.0) - prev) * (Type(1.0) - specindex) * (Type(1.0) - specref);
        Type p10 = prev * (sensindex * (Type(1.0) - sensref)) + (Type(1.0) - prev) * (Type(1.0) - specindex) * specref;
        Type p01 = prev * (Type(1.0) - sensindex) * sensref   + (Type(1.0) - prev) * specindex * (Type(1.0) - specref);
        Type p00 = prev * (Type(1.0) - sensindex) * (Type(1.0) - sensref) + (Type(1.0) - prev) * (specindex * specref);

        p11 = CppAD::CondExpLt(p11, eps, eps, p11);
        p10 = CppAD::CondExpLt(p10, eps, eps, p10);
        p01 = CppAD::CondExpLt(p01, eps, eps, p01);
        p00 = CppAD::CondExpLt(p00, eps, eps, p00);

        /* ===== LIKELIHOOD CONTTRIBUTIONS ===== */
        nll -= y11(i) * log(p11) + y10(i) * log(p10) + y01(i) * log(p01) + y00(i) * log(p00);

    }

    /* =====  SROC calculations ===== */
    vector<Type> logitsens(nspec);
    vector<Type> sens(nspec);
    Type exp_beta      = exp(-beta);
    Type exp_beta_half = exp(-beta / Type(2.0));
    
    for(int i = 0; i < nspec; i++)
    {
        Type sp = spec(i);
        logitsens(i) = Lambda * exp_beta_half - exp_beta * logit(sp);
        sens(i) = invlogit(logitsens(i));
    }

    /* =====  REPORTS ===== */
    Type sigma2_alpha = sigma_alpha * sigma_alpha;
    Type sigma2_theta = sigma_theta * sigma_theta;
    Type sigma2_prev  = sigma_prev  * sigma_prev;


    REPORT(mu_prev);
    REPORT(Lambda);
    REPORT(Theta);
    REPORT(beta);
    REPORT(sigma2_prev);
    REPORT(sigma2_alpha);
    REPORT(sigma2_theta);
    REPORT(mu_A_ref);
    REPORT(mu_B_ref);
    REPORT(logitsens);
    REPORT(sens);
    
    ADREPORT(mu_prev);
    ADREPORT(Lambda);
    ADREPORT(Theta);
    ADREPORT(beta);
    ADREPORT(sigma2_prev);
    ADREPORT(sigma2_alpha);
    ADREPORT(sigma2_theta);
    ADREPORT(mu_A_ref);
    ADREPORT(mu_B_ref);
    ADREPORT(logitsens);
    ADREPORT(sens);

    return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this