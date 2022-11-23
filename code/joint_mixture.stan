functions {
    /* GAUSSIAN @rewritten from Pfeuffer 2017 / Pinkney 2021 @param rho (-1, 1) */
    real normal_copula(real u, real v, real rho) {
        real rho_sq = square(rho);
        real s = inv_Phi(u);
        real t = inv_Phi(v);
        return -0.5 * log1m(rho_sq) - ((rho_sq * (square(s) + square(t)) - (2. * rho * s * t)) / (2. * (1. - rho_sq)));
    }

    /* GUMBEL @written by Ben Goodrich, 2017 (https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ) @param theta */
    real gumbel_copula(real u, real v, real theta) {
        real neg_log_u = -log(u);
        real log_neg_log_u = log(neg_log_u);
        real neg_log_v = -log(v);
        real log_neg_log_v = log(neg_log_v);
        real log_temp = log_sum_exp(theta * log_neg_log_u, theta * log_neg_log_v);
        real theta_m1 = theta - 1;

        if (theta < 1) reject("theta must be >= 1");
        if (is_inf(theta))
            return (u == v) ? 0 : negative_infinity();

        return theta_m1 * log_neg_log_u + theta_m1 * log_neg_log_v + neg_log_u + neg_log_v - exp(log_temp / theta) +
              log_sum_exp(2. * theta_m1  / -theta * log_temp, log(theta_m1) + (1. - 2. * theta) / theta * log_temp);
     }

    /* CLAYTON @written by Andre Pfeuffer (https://groups.google.com/g/stan-users/c/hnUtkMYlLhQ/m/XdX3u1vDAAAJ) @param theta */
    real clayton_copula(real u, real v, real theta) {
        return log1p(theta) - (theta + 1) * (log(u) + log(v)) - (1 + 2 * theta) / theta * log(pow(u, -theta) + pow(v, -theta) - 1);
    }
}

data {
    int<lower=0> N; // number of records
    vector[N] EL;
    vector<lower=EL>[N] ER;
    vector[N] CL;
    vector<lower=CL>[N] CR;
    vector[N] S1L;
    vector<lower=S1L>[N] S1R;

    vector<lower=0>[2] mu_mean_prior;
    vector<lower=0>[2] mu_sigma_prior;

    vector<lower=0>[2] par1_mean_prior;
    vector<lower=0>[2] par1_sigma_prior;
}

transformed data {
    int D = 3; // number of distributions for each GI and IP (Gamma, Weibull, Lognormal)
    int K = 4; // number of copulas (Independence, Gaussian, Gumbel, Clayton)
    int M = D * D * K; // total number of combinations
}

parameters {
    vector<lower=0, upper=1>[N] e_raw; // uniform normalized value for sampling between EL and ER
    vector<lower=0, upper=1>[N] c_raw; // uniform normalized value for sampling between CL and CR
    vector<lower=0, upper=1>[N] s1_raw; // uniform normalized value for sampling between S1L and S1R

    real theta_normal_raw;

    vector<lower=0>[2] mu; //mean of GI and IP
    vector<lower=0>[2] param1_weibull;

    simplex[M] lambda; // mixing proportions
}

transformed parameters {
    vector<lower = 0>[2] sigma; // SDs of GI and IP

    real theta_normal = inv_logit(theta_normal_raw); // theta_normal is between zero and one
    //real tau = asin(2. * theta_normal / pi()); // Kendalls tau
    real tau = 2. / pi() * asin(theta_normal); // from this formula

    vector[D] param1[2];
    vector[D] param2[2];

    vector[M] lps = log(lambda); // likelihood
    {
        real theta_clayton = 2. * tau / (1. - tau);
        real theta_gumbel = 1. / (1. - tau);
        real theta_independence = 1.0;

        vector[N] e; // infection of infector
        vector[N] c; // infection of infectee
        vector[N] s1; // symptoms of infector

        e = EL + (ER - EL) .* e_raw;
        for (n in 1:N) {
            c[n] = (CL[n] < e[n]) ? fma(CR[n] - e[n], c_raw[n], e[n]) : fma(CR[n] - CL[n], c_raw[n], CL[n]);
            s1[n] = (S1L[n] < e[n]) ? fma(S1R[n] - e[n], s1_raw[n], e[n]) : fma(S1R[n] - S1L[n], s1_raw[n], S1L[n]);
        }

        for (j in 1:2) {  // 1 is GI, 2 is IP
            // Weibull distribution
            param1[j, 2] = param1_weibull[j];
            param2[j, 2] = mu[j] / tgamma(1.0 + 1.0 / param1[j, 2]);
            sigma[j] = param2[j, 2] * sqrt(tgamma(1.0 + 2.0 / param1[j, 2]) - square(tgamma(1.0 + 1.0 / param1[j, 2])));

            // Gamma distribution
            param1[j, 1] = square(mu[j] / sigma[j]);
            param2[j, 1] = mu[j] / square(sigma[j]);

            // Lognormal distribution
            param2[j, 3] = sqrt(log(square(sigma[j] / mu[j]) + 1.0));
            param1[j, 3] = log(mu[j]) - square(param2[j, 3]) / 2.0;
        }

        int i; int j; int k;
        for (n in 1:M) {
            // transforming 1d array to 3d array (n -> [i, j, k])
            // see https://stackoverflow.com/questions/11316490/convert-a-1d-array-index-to-a-3d-array-index
            k = (n - 1) % K + 1;
            j = ((n - 1) / K) % D + 1;
            i = (n - 1) / (D * K) + 1;

            vector[N] x; vector[N] y; // copula vectors
            // Generation interval contribution
            if (i == 1) {
                lps[n] += gamma_lpdf(c - e | param1[1, i], param2[1, i]);
                for (l in 1:N)
                    x[l] = gamma_cdf(c[l] - e[l], param1[1, i], param2[1, i]);
            } else if (i == 2) {
                lps[n] += weibull_lpdf(c - e | param1[1, i], param2[1, i]);
                for (l in 1:N)
                    x[l] = weibull_cdf(c[l] - e[l], param1[1, i], param2[1, i]);
            } else {
                lps[n] += lognormal_lpdf(c - e | param1[1, i], param2[1, i]);
                for (l in 1:N)
                    x[l] = lognormal_cdf(c[l] - e[l], param1[1, i], param2[1, i]);
            }
            // Incubation period contribution
            if (j == 1) {
                lps[n] += gamma_lpdf(s1 - e | param1[2, j], param2[2, j]);
                for (l in 1:N)
                    y[l] = gamma_cdf(s1[l] - e[l], param1[2, j], param2[2, j]);
            } else if (j == 2) {
                lps[n] += weibull_lpdf(s1 - e | param1[2, j], param2[2, j]);
                for (l in 1:N)
                    y[l] = weibull_cdf(s1[l] - e[l], param1[2, j], param2[2, j]);
            } else {
                lps[n] += lognormal_lpdf(s1 - e | param1[2, j], param2[2, j]);
                for (l in 1:N)
                    y[l] = lognormal_cdf(s1[l] - e[l], param1[2, j], param2[2, j]);
            }
            // Copula contribution
            if (k == 1)
                for (l in 1:N)
                    lps[n] += gumbel_copula(x[l], y[l], theta_independence);
            else if (k == 2)
                for (l in 1:N)
                    lps[n] += normal_copula(x[l], y[l], theta_normal);
            else if (k == 3)
                for (l in 1:N)
                    lps[n] += gumbel_copula(x[l], y[l], theta_gumbel);
            else
                for (l in 1:N)
                    lps[n] += clayton_copula(x[l], y[l], theta_clayton);
        }
    }
}

model {
    /* priors */
    mu ~ normal(mu_mean_prior, mu_sigma_prior);
    param1_weibull ~ normal(par1_mean_prior, par1_sigma_prior);

    e_raw ~ normal(.5, 2.0);
    c_raw ~ normal(.5, 2.0);
    s1_raw ~ normal(.5, 2.0);

    theta_normal_raw ~ logistic(0.0, 1.0);

    target += log_sum_exp(lps);
}

generated quantities {
    matrix<lower = 0, upper = 1>[D, K] q[D];
    int index_GI; int index_IP; int index_copula;
    {
        vector[D * D * K] q_ = exp(lps - log_sum_exp(lps)); // see comments at https://statmodeling.stat.columbia.edu/2017/08/21/mixture-models-stan-can-use-log_mix/
        int i; int j; int k;
        for (n in 1:(D * D * K)) {
            k = (n - 1) % K + 1;
            j = ((n - 1) / K) % D + 1;
            i = (n - 1) / (D * K) + 1;
            q[i, j, k] = q_[n];
        }
        int n_selected = categorical_rng(q_);
        index_copula = (n_selected - 1) % K + 1;
        index_IP = ((n_selected - 1) / K) % D + 1;
        index_GI = (n_selected - 1) / (D * K) + 1;
    }
}
