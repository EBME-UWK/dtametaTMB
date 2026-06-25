#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type RutterGatsonis(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y);
    DATA_VECTOR(n);
    DATA_VECTOR(x);
    DATA_FACTOR(study);
    DATA_VECTOR(spec);

    int N     = y.size();
    int nspec = spec.size();
    
    /* ===== VALIDATE x ===== */
    for(int i = 0; i < x.size(); i++)
    {
        Type xi = x(i);

        if(!R_finite(asDouble(xi)))
           error("x contains non-finite values");

        if(!(xi == Type(-0.5) || xi == Type(0.5)))
           error("x must contain only -0.5 or 0.5");
    }

    /* =====  VALIDATE SPEC GRID ===== */
    for(int i=0; i<nspec; i++)
    {
        double sp = asDouble(spec(i));
        if(!R_finite(sp)) error("spec contains non-finite values");
        if(sp <= 0.0 || sp >= 1.0) error("All values in spec must satisfy 0 < spec < 1");
    }

    /* ===== FIXED EFFECTS ===== */
    PARAMETER(Lambda);
    PARAMETER(Theta);
    PARAMETER(beta);
    PARAMETER(log_sigma_alpha);
    PARAMETER(log_sigma_theta);

    Type sigma_alpha = exp(log_sigma_alpha);
    Type sigma_theta = exp(log_sigma_theta);

    /* =====  RANDOM EFFECTS (per study) ===== */
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(theta);
    
    /* =====  NEGATIVE LOG-LIKELIHOOD ===== */
    Type nll = 0.0;

    /* =====  Random effects (vectorized) ===== */
    nll -= sum(dnorm(alpha, Type(0.0), sigma_alpha, true));
    nll -= sum(dnorm(theta, Type(0.0), sigma_theta, true));

    /* =====  Binomial likelihood ===== */
    for(int i = 0; i < N; i++)
    {
        int s = study(i);   //
        Type eta = ( Theta + theta(s) + ( Lambda + alpha(s) ) * x(i) ) * exp( -beta * x(i) );
        Type p = invlogit(eta);
        nll -= dbinom(y(i), n(i), p, true);
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

    /* =====  REPORT (R output) ===== */
    Type sigma2_alpha = sigma_alpha * sigma_alpha;
    Type sigma2_theta = sigma_theta * sigma_theta;

    REPORT(Lambda);
    REPORT(Theta);
    REPORT(beta);
    REPORT(sigma2_alpha);
    REPORT(sigma2_theta);
    REPORT(logitsens);
    REPORT(sens);
  
    ADREPORT(Lambda);
    ADREPORT(Theta);
    ADREPORT(beta);
    ADREPORT(sigma2_alpha);
    ADREPORT(sigma2_theta);
    ADREPORT(logitsens);
    ADREPORT(sens);

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this