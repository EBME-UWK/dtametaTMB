#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type RutterGatsonisReg(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */
    DATA_VECTOR(y);
    DATA_VECTOR(n);
    DATA_VECTOR(x);
    DATA_MATRIX(Z);      // covariate matrix
    DATA_MATRIX(Z_pred); // prediction matrix
    DATA_FACTOR(study);
    DATA_VECTOR(spec);

    int N      = y.size();
    int nspec  = spec.size();
    int ngroup = Z_pred.rows();
    int p = Z.cols();

    /* ===== INPUT CHECKS ===== */

    if(n.size() != N)
       error("n and y must have same length");
    if(x.size() != N)
       error("x and y must have same length");
    if(study.size() != N)
       error("study and y must have same length");
    if(Z.rows() != N)
       error("Z must have N rows");
    if(Z_pred.cols() != Z.cols())
       error("Z_pred and Z must have the same number of columns");
    if(Z_pred.rows() < 1 || Z_pred.cols() < 1)
        error("Z_pred must contain at least one row and one column");

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
    PARAMETER_VECTOR(accuracy_coef);
    PARAMETER_VECTOR(threshold_coef);
    PARAMETER_VECTOR(shape_coef);

    if(accuracy_coef.size() != p)
       error("Length of accuracy_coef must equal ncol(Z)");
    if(threshold_coef.size() != p)
       error("Length of threshold_coef must equal ncol(Z)");
    if(shape_coef.size() != p)
       error("Length of shape_coef must equal ncol(Z)");

    vector<Type> accuracy_vec  = Z * accuracy_coef;
    vector<Type> threshold_vec = Z * threshold_coef;   
    vector<Type> shape_vec     = Z * shape_coef;

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
        Type eta = ( threshold_vec(i) + theta(s) + ( accuracy_vec(i) + alpha(s) ) * x(i) ) * exp( -shape_vec(i) * x(i) );
        Type p = invlogit(eta);
        nll -= dbinom(y(i), n(i), p, true);
    }

    /* =====  SROC calculations ===== */
    matrix<Type> logitsens(nspec,ngroup);
    matrix<Type> sens(nspec,ngroup);

    vector<Type> Lambda_Pred = Z_pred * accuracy_coef;
    vector<Type> Theta_Pred  = Z_pred * threshold_coef;
    vector<Type> beta_Pred   = Z_pred * shape_coef;

    for(int g=0; g<ngroup; g++)
    {
        for(int i=0; i<nspec; i++)
        {
            logitsens(i,g) = Lambda_Pred(g) * exp(-beta_Pred(g) / Type(2.0)) - exp(-beta_Pred(g)) * logit(spec(i));
            sens(i,g) = invlogit(logitsens(i,g));
        }
    }

    /* =====  REPORT (R output) ===== */
    Type sigma2_alpha = sigma_alpha * sigma_alpha;
    Type sigma2_theta = sigma_theta * sigma_theta;

    REPORT(accuracy_coef);
    REPORT(threshold_coef);
    REPORT(shape_coef);
    REPORT(sigma2_alpha);
    REPORT(sigma2_theta);
    REPORT(Lambda_Pred);
    REPORT(Theta_Pred);
    REPORT(beta_Pred);
    REPORT(logitsens);
    REPORT(sens);
  
    ADREPORT(accuracy_coef);
    ADREPORT(threshold_coef);
    ADREPORT(shape_coef);
    ADREPORT(sigma2_alpha);
    ADREPORT(sigma2_theta);
    ADREPORT(Lambda_Pred);
    ADREPORT(Theta_Pred);
    ADREPORT(beta_Pred);
    ADREPORT(logitsens);
    ADREPORT(sens);

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this