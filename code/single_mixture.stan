data {
    int<lower=0> N; // number of records
    vector[N] EL;
    vector<lower=EL>[N] ER;
    vector[N] CL;
    vector<lower=CL>[N] CR;

    real<lower=0> mu_mean_prior;
    real<lower=0> mu_sigma_prior;

    real<lower=0> par1_mean_prior;
    real<lower=0> par1_sigma_prior;
}

transformed data {
    int D = 3; // number of distributions (Gamma, Weibull, Lognormal)
}

parameters {
    vector<lower=0, upper=1>[N] e_raw; // uniform normalized value for sampling between EL and ER
    vector<lower=0, upper=1>[N] c_raw; // uniform normalized value for sampling between CL and CR

    real<lower=0> mu; //mean of distribution
    real<lower=0> param1_weibull;

    simplex[D] weight; // mixing proportions
}

transformed parameters {
    real<lower=0> sigma; // SD of distribution

    vector[D] param1; // distribution parameters
    vector[D] param2; // distribution parameters

    vector[D] lps = log(weight); // likelihood
    {
        vector[N] e; // infection of infector
        vector[N] c; // infection of infectee

        e = EL + (ER - EL) .* e_raw;
        for (n in 1:N) {
            c[n] = (CL[n] < e[n]) ? fma(CR[n] - e[n], c_raw[n], e[n]) : fma(CR[n] - CL[n], c_raw[n], CL[n]);
        }

            // Weibull distribution
            param1[2] = param1_weibull;
            param2[2] = mu / tgamma(1.0 + 1.0 / param1[2]);
            sigma = param2[2] * sqrt(tgamma(1.0 + 2.0 / param1[2]) - square(tgamma(1.0 + 1.0 / param1[2])));

            // Gamma distribution
            param1[1] = square(mu / sigma);
            param2[1] = mu / square(sigma);

            // Lognormal distribution
            param2[3] = sqrt(log(square(sigma / mu) + 1.0));
            param1[3] = log(mu) - square(param2[3]) / 2.0;

        int i;
        for (n in 1:D) {
            i = (n - 1) % D + 1;

            if (i == 1) {
                lps[n] += gamma_lpdf(c - e | param1[i], param2[i]);
            } else if (i == 2) {
                lps[n] += weibull_lpdf(c - e | param1[i], param2[i]);
            } else {
                lps[n] += lognormal_lpdf(c - e | param1[i], param2[i]);
            }
        }
    }
}

model {
    /* priors */
    mu ~ normal(mu_mean_prior, mu_sigma_prior);
    param1_weibull ~ normal(par1_mean_prior, par1_sigma_prior);

    e_raw ~ normal(.5, 2.0);
    c_raw ~ normal(.5, 2.0);

    target += log_sum_exp(lps);
}

generated quantities {
    vector<lower = 0, upper = 1>[D] q;
    int index_dist;
    {
        vector[D] q_ = exp(lps - log_sum_exp(lps)); // see comments at https://statmodeling.stat.columbia.edu/2017/08/21/mixture-models-stan-can-use-log_mix/
        int i;
        for (n in 1:D) {
            i = (n - 1) % D + 1;
            q[i] = q_[n];
        }
        int n_selected = categorical_rng(q_);
        index_dist = (n_selected - 1) % D + 1;
    }
}
