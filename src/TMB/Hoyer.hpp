#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* ===== Log survival function ===== */

template<class Type>
Type safe_logspace_sub(Type a, 
                       Type b) {
    Type eps = Type(1e-12);
    return CppAD::CondExpLt(a - b, eps, a + log(eps), logspace_sub(a, b));
}

template<class Type>
Type safe_time(Type t) {
    // Prevent log(0) or log(negative)
    return CppAD::CondExpGt(t, Type(1e-12), t, Type(1e-12));
}

template<class Type>
Type log_survival(Type t,
                  Type eta,
                  Type lambda,
                  int dist){

    t = safe_time(t);

    if(dist == 1) {            // Weibull
        Type log_a = (log(t) - eta) / lambda;
        Type log_a_clamped = CppAD::CondExpGt(log_a, Type(50), Type(50), log_a);
        Type a = exp(log_a_clamped);
        return -a;

    } else if(dist == 2) {     // Log-normal
        Type z = (log(t) - eta) / lambda;
        Type Phi = pnorm(z);
        Type eps = Type(1e-16);
        Phi = CppAD::CondExpLt(Phi, eps, eps, Phi);
        Type logPhi = log(Phi);
        return safe_logspace_sub(Type(0.0), logPhi);

    } else if(dist == 3) {     // Log-logistic
        Type z = (eta - log(t)) / lambda;
        return -logspace_add(Type(0.0), -z);

    } else {
        error("dist must be 1 (=weibull), 2 (=lognormal) or 3 (=loglogistic)");
        return Type(0);
    }
}

/* ===== logit( P(T > t) ) computed in log-space ===== */

template<class Type>
Type logit_survival(Type t,
                    Type eta,
                    Type lambda,
                    int dist){

    // --- Ensure valid time ---
    t = safe_time(t);

    if(dist == 1) {            // Weibull
        Type log_a = (log(t) - eta) / lambda;
        Type log_a_clamped = CppAD::CondExpGt(log_a, Type(50), Type(50), log_a);
        Type a = exp(log_a_clamped);
        Type logS = -a;
        Type log1mS = safe_logspace_sub(Type(0.0), logS);
        return logS - log1mS;

    } else if(dist == 2) {     // Log-normal
        Type z = (log(t) - eta) / lambda;
        Type Phi = pnorm(z);
        Type eps = Type(1e-16);
        Phi = CppAD::CondExpLt(Phi, eps, eps, Phi);
        Type logPhi = log(Phi);
        Type logS = safe_logspace_sub(Type(0.0), logPhi);
        return logS - logPhi;

    } else if(dist == 3) {     // Log-logistic
        return (eta - log(t)) / lambda;

    } else {
        error("dist must be 1 (=weibull), 2 (=lognormal) or 3 (=loglogistic)");
        return Type(0);
    }
}

/* ===== Main objective ===== */

template<class Type>
Type Hoyer(objective_function<Type>* obj)
{
    using namespace density;

    /* ===== DATA ===== */

    DATA_VECTOR(lowerB);
    DATA_VECTOR(upperB);
    DATA_VECTOR(events0);
    DATA_VECTOR(events1);
    DATA_VECTOR(threshold);
    DATA_IVECTOR(ctype);
    DATA_IVECTOR(study);
    DATA_INTEGER(nstudy);
    DATA_INTEGER(dist);

    int n = lowerB.size();

    if(upperB.size()  != n) error("upperB size mismatch");
    if(events0.size() != n) error("events0 size mismatch");
    if(events1.size() != n) error("events1 size mismatch");
    if(ctype.size()   != n) error("ctype size mismatch");
    if(study.size()   != n) error("study size mismatch");

    /* ===== FIXED EFFECTS ===== */

    PARAMETER(beta0);
    PARAMETER(log_lambda0);
    PARAMETER(beta1);
    PARAMETER(log_lambda1);

    Type lambda0 = exp(log_lambda0);
    Type lambda1 = exp(log_lambda1);

    /* ===== RANDOM EFFECTS ===== */

    PARAMETER(log_su0);
    PARAMETER(log_su1);
    PARAMETER(rho_trans);

    Type su0 = exp(log_su0);
    Type su1 = exp(log_su1);

    // slight shrink to avoid exact singularity
    Type rho = Type(0.9999) * tanh(rho_trans);

    Type covu0u1 = rho * su0 * su1;

    PARAMETER_VECTOR(u0);
    PARAMETER_VECTOR(u1);

    if(u0.size() != nstudy || u1.size() != nstudy)
        error("Random effect length mismatch");

    matrix<Type> Sigma(2,2);
    Sigma(0,0) = su0 * su0;
    Sigma(1,1) = su1 * su1;
    Sigma(0,1) = covu0u1;
    Sigma(1,0) = covu0u1;

    MVNORM_t<Type> neg_log_density(Sigma);

    /* ===== NEGATIVE LOG-LIKELIHOOD ===== */

    Type nll = 0;

    /* ===== RANDOM EFFECT CONTRIBUTION ===== */

    vector<Type> uj(2);

    for(int j=0; j<nstudy; j++) {
        uj(0) = u0(j);
        uj(1) = u1(j);
        nll += neg_log_density(uj);
    }

    /* ===== PRECOMPUTE STUDY ETAS ===== */

    vector<Type> eta0_study(nstudy);
    vector<Type> eta1_study(nstudy);

    for(int j=0; j<nstudy; j++) {
        eta0_study(j) = beta0 + u0(j);
        eta1_study(j) = beta1 + u1(j);
    }

    /* ===== OBSERVATION CONTRIBUTION ===== */

    for(int i=0; i<n; i++) {

        int s = study(i);
        if(s < 0 || s >= nstudy) error("study index out of bounds");

        Type eta0 = eta0_study(s);
        Type eta1 = eta1_study(s);

        Type logp0;
        Type logp1;

        /* ===== LEFT CENSORED ===== */

        if(ctype(i) == 1) {
            Type logS0u = log_survival(upperB(i), eta0, lambda0, dist);
            Type logS1u = log_survival(upperB(i), eta1, lambda1, dist);
            logp0 = safe_logspace_sub(Type(0.0), logS0u);
            logp1 = safe_logspace_sub(Type(0.0), logS1u);

        /* ===== INTERVAL CENSORED ===== */

        } else if(ctype(i) == 2) {
            Type logS0l = log_survival(lowerB(i), eta0, lambda0, dist);
            Type logS0u = log_survival(upperB(i), eta0, lambda0, dist);
            Type logS1l = log_survival(lowerB(i), eta1, lambda1, dist);
            Type logS1u = log_survival(upperB(i), eta1, lambda1, dist);
            logp0 = safe_logspace_sub(logS0l, logS0u);
            logp1 = safe_logspace_sub(logS1l, logS1u);

        /* ===== RIGHT CENSORED ===== */

        } else if(ctype(i) == 3) {
            logp0 = log_survival(lowerB(i), eta0, lambda0, dist);
            logp1 = log_survival(lowerB(i), eta1, lambda1, dist);

        } else {
            error("ctype must be 1, 2 or 3");
        }

        nll -= events0(i) * logp0;
        nll -= events1(i) * logp1;
    }

    /* ===== SUMMARY ROC ===== */

    const int nt = threshold.size();

    vector<Type> logitSurv1(nt);
    vector<Type> logitSurv0(nt);

    for(int k=0; k<nt; k++) {
        logitSurv1(k) = logit_survival(threshold(k), beta1, lambda1, dist);
        logitSurv0(k) = logit_survival(threshold(k), beta0, lambda0, dist);
    }

    /* ===== FINAL SAFETY GUARD ===== */

    if(!CppAD::isfinite(nll)) return Type(1e20);

    /* ===== REPORT ===== */

    REPORT(su0);
    REPORT(su1);
    REPORT(rho);
    REPORT(covu0u1);

    REPORT(logitSurv1); //=logitSens for testdirection="greater" but logitfnr  for testdirection="less"
    REPORT(logitSurv0); //=logitfpr  for testdirection="greater" but logitspec for testdirection="less"

    ADREPORT(beta0);
    ADREPORT(lambda0);
    ADREPORT(su0);
    ADREPORT(beta1);
    ADREPORT(lambda1);
    ADREPORT(su1);
    ADREPORT(rho);
    ADREPORT(covu0u1);
    ADREPORT(logitSurv1); //=logitSens for testdirection="greater" but logitfnr  for testdirection="less"
    ADREPORT(logitSurv0); //=logitfpr  for testdirection="greater" but logitspec for testdirection="less"

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this