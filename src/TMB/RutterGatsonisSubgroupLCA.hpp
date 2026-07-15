#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type RutterGatsonisSubgroupLCA(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y11);
    DATA_VECTOR(y10);
    DATA_VECTOR(y01);
    DATA_VECTOR(y00);
    DATA_MATRIX(Z);
    DATA_MATRIX(Z_pred);
    DATA_VECTOR(spec);

    int N      = y11.size();
    int nspec  = spec.size();
    int ngroup = Z_pred.rows();
    int p      = Z.cols();

    if(y10.size() != N)
        error("y10 and y11 must have equal lengths");

    if(y01.size() != N)
        error("y01 and y11 must have equal lengths");

    if(y00.size() != N)
        error("y00 and y11 must have equal lengths");

    if(Z.rows() != N)
       error("Z must equal N rows");
    if(Z_pred.cols() != Z.cols())
       error("Z_pred and Z must have the same number of columns");
    if(Z_pred.rows() < 1 || Z_pred.cols() < 1)
        error("Z_pred must contain at least one row and one column");
        
    /* =====  VALIDATE SPEC GRID ===== */
    for(int i=0; i<nspec; i++)
    {
        double sp = asDouble(spec(i));
        if(!R_finite(sp)) error("spec contains non-finite values");
        if(sp <= 0.0 || sp >= 1.0) error("All values in spec must satisfy 0 < spec < 1");
    }
    
    /* ===== FIXED EFFECTS ===== */
    PARAMETER_VECTOR(prev_coef); // logit prevalence
    PARAMETER_VECTOR(accuracy_coef);
    PARAMETER_VECTOR(threshold_coef);
    PARAMETER_VECTOR(shape_coef);

    if(prev_coef.size() != p)
       error("Length of prev_coef must equal ncol(Z)");
    if(accuracy_coef.size() != p)
       error("Length of accuracy_coef must equal ncol(Z)");
    if(threshold_coef.size() != p)
       error("Length of threshold_coef must equal ncol(Z)");
    if(shape_coef.size() != p)
       error("Length of shape_coef must equal ncol(Z)");

    vector<Type> prev_vec      = Z * prev_coef;
    vector<Type> accuracy_vec  = Z * accuracy_coef;
    vector<Type> threshold_vec = Z * threshold_coef;
    vector<Type> shape_vec     = Z * shape_coef;

    PARAMETER_VECTOR(log_sigma_prev_coef);

    if(log_sigma_prev_coef.size() != p)
        error("log_sigma_prev_coef length must equal number ncol(Z)");

    vector<Type> log_sigma_prev_vec = Z * log_sigma_prev_coef;
    vector<Type> sigma_prev_vec(N);
    for(int i = 0; i < N; i++)
    {
       sigma_prev_vec(i) = exp(log_sigma_prev_vec(i));
    }

    PARAMETER(log_sigma_alpha);
    PARAMETER(log_sigma_theta);

    Type sigma_alpha  = exp(log_sigma_alpha);
    Type sigma_theta  = exp(log_sigma_theta);
    
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
    for(int i = 0; i < N; i++)
    {
        nll -= dnorm(prevu(i), Type(0.0), sigma_prev_vec(i), true);
    }

    for(int i = 0; i < N; i++)
    {   
        /* ===== LINEAR PREDICTORS AND TRANSFORM ===== */
        Type prev = invlogit( prev_vec(i) + prevu(i) );
        Type etaD = ( threshold_vec(i) + theta(i) + ( accuracy_vec(i) + alpha(i) )/Type(2.0) ) * exp( -shape_vec(i)/Type(2.0) );
        Type etaH = ( threshold_vec(i) + theta(i) - ( accuracy_vec(i) + alpha(i) )/Type(2.0) ) * exp(  shape_vec(i)/Type(2.0) );
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

        /* ===== LIKELIHOOD CONTRIBUTIONS ===== */
        nll -= y11(i) * log(p11) + y10(i) * log(p10) + y01(i) * log(p01) + y00(i) * log(p00);

    }

 /* =====  SROC calculations ===== */
    matrix<Type> logitsens(nspec,ngroup);
    matrix<Type> sens(nspec,ngroup);

    vector<Type> Lambda_Pred = Z_pred * accuracy_coef;
    vector<Type> Theta_Pred  = Z_pred * threshold_coef;
    vector<Type> beta_Pred   = Z_pred * shape_coef;
    vector<Type> prev_Pred   = Z_pred * prev_coef;

    for(int g=0; g<ngroup; g++)
    {
        for(int i=0; i<nspec; i++)
        {
            logitsens(i,g) = Lambda_Pred(g) * exp(-beta_Pred(g) / Type(2.0)) - exp(-beta_Pred(g)) * logit(spec(i));
            sens(i,g) = invlogit(logitsens(i,g));
        }
    }

    /* =====  REPORTS ===== */
    Type sigma2_alpha = sigma_alpha * sigma_alpha;
    Type sigma2_theta = sigma_theta * sigma_theta;
    vector<Type> sigma2_prev(ngroup);
    vector<Type> log_sigma_prev_pred = Z_pred * log_sigma_prev_coef;
    for(int g=0; g<ngroup; g++)
    {  
        sigma2_prev(g) = exp(Type(2.0) * log_sigma_prev_pred(g));
    }

    REPORT(prev_coef);
    REPORT(accuracy_coef);
    REPORT(threshold_coef);
    REPORT(shape_coef);
    REPORT(sigma2_prev);
    REPORT(sigma2_alpha);
    REPORT(sigma2_theta);
    REPORT(Lambda_Pred);
    REPORT(Theta_Pred);
    REPORT(beta_Pred);
    REPORT(prev_Pred);
    REPORT(mu_A_ref);
    REPORT(mu_B_ref);
    REPORT(logitsens);
    REPORT(sens);
  
    ADREPORT(prev_coef);
    ADREPORT(accuracy_coef);
    ADREPORT(threshold_coef);
    ADREPORT(shape_coef);
    ADREPORT(sigma2_prev);
    ADREPORT(sigma2_alpha);
    ADREPORT(sigma2_theta);
    ADREPORT(Lambda_Pred);
    ADREPORT(Theta_Pred);
    ADREPORT(beta_Pred);
    ADREPORT(prev_Pred); 
    ADREPORT(mu_A_ref);
    ADREPORT(mu_B_ref);
    ADREPORT(logitsens);
    ADREPORT(sens);

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this