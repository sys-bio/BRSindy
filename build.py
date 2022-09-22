
# creation of ansatz reactions

from itertools import combinations, permutations, product


def model(n_species, rxn_types=None):

    species = []
    spec = []
    for i in range(n_species):
        species.append('S' + str(i))
        spec.append(i)

    if rxn_types is None:
        rxn_types = ['syn', 'deg', 'deg2', 'uni-uni', 'bi-uni', 'uni-bi', 'bi-bi']

    bi_species = None
    bi_spec = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types or 'bi-bi' in rxn_types:
        bi_species = list(combinations(species, 2))
        bi_spec = list(combinations(spec, 2))
        bi_species_ident = []
        bi_spec_ident = []
        for each in species:
            bi_species_ident.append((each, each))
        for each in spec:
            bi_spec_ident.append((each, each))
        bi_species = bi_species_ident + bi_species
        bi_spec = bi_spec_ident + bi_spec

    syn_r = []
    deg_r = []
    deg2_r = []
    uni_uni_r = []
    bi_uni_r = []
    uni_bi_r = []
    bi_bi_r = []

    function_count = 0
    mat = [[] for _ in range(n_species)]
    kin = []

    if 'syn' in rxn_types:
        for i, each in enumerate(spec):
            kin.append([])
            function_count += 1
            for j, item in enumerate(mat):
                if each == j:
                    mat[j].append(1)
                else:
                    mat[j].append(0)
        for i, each in enumerate(species):
            syn_r.append('-> ' + each + ';')

    if 'deg' in rxn_types:
        for i, each in enumerate(spec):
            kin.append([])
            function_count += 1
            for j, item in enumerate(mat):
                if each == j:
                    mat[j].append(-1)
                    kin[-1].append(each)
                else:
                    mat[j].append(0)
        for i, each in enumerate(species):
            deg_r.append(each + ' -> ;')

    if 'deg2' in rxn_types:
        for i, each in enumerate(bi_spec):
            kin.append([])
            function_count += 1
            for j, item in enumerate(mat):
                mat[j].append(0)
                for k, every in enumerate(bi_spec[i]):
                    if every == j:
                        mat[j][-1] -= 1
                        kin[-1].append(every)
        for each in bi_species:
            deg2_r.append(each[0] + ' + ' + each[1] + ' -> ;')

    if 'uni-uni' in rxn_types:
        uni_permutations = list(permutations(species, 2))
        uni_perm = list(permutations(spec, 2))
        for i, each in enumerate(uni_perm):
            uni_perm[i] = list(uni_perm[i])
            kin.append([uni_perm[i][0]])
            function_count += 1
            for j, item in enumerate(mat):
                if uni_perm[i][0] == j:
                    mat[j].append(-1)
                if uni_perm[i][1] == j:
                    mat[j].append(1)
        for each in uni_permutations:
            uni_uni_r.append(each[0] + ' -> ' + each[1] + ';')

    tri_species = None
    tri_spec = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types:
        tri_species = list(product(species, bi_species))
        tri_spec = list(product(spec, bi_spec))
        for i, each in enumerate(tri_spec):
            tri_spec[i] = list(tri_spec[i])
            tri_spec[i][1] = list(tri_spec[i][1])

    if 'bi-uni' in rxn_types:
        for i, each in enumerate(tri_spec):
            function_count += 1
            kin.append([])
            kin[-1].append(tri_spec[i][1][0])
            kin[-1].append(tri_spec[i][1][1])
            for j, item in enumerate(mat):
                mat[j].append(0)
                if tri_spec[i][0] == j:
                    mat[j][-1] += 1
                if tri_spec[i][1][0] == j:
                    mat[j][-1] -= 1
                if tri_spec[i][1][1] == j:
                    mat[j][-1] -= 1

        for each in tri_species:
            bi_uni_r.append(each[1][0] + ' + ' + each[1][1] + ' -> ' + each[0] + ';')

    if 'uni-bi' in rxn_types:
        for i, each in enumerate(tri_spec):
            function_count += 1
            kin.append([])
            kin[-1].append(tri_spec[i][0])
            for j, item in enumerate(mat):
                mat[j].append(0)
                if tri_spec[i][0] == j:
                    mat[j][-1] -= 1
                if tri_spec[i][1][0] == j:
                    mat[j][-1] += 1
                if tri_spec[i][1][1] == j:
                    mat[j][-1] += 1

        for each in tri_species:
            uni_bi_r.append(each[0] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')

    if 'bi-bi' in rxn_types:
        quad_species = list(permutations(bi_species, 2))
        quad_spec = list(permutations(bi_spec, 2))
        for i, each in enumerate(quad_spec):
            quad_spec[i] = list(quad_spec[i])
            quad_spec[i][0] = list(quad_spec[i][0])
            quad_spec[i][1] = list(quad_spec[i][1])
        for i, each in enumerate(quad_spec):
            function_count += 1
            kin.append([])
            kin[-1].append(quad_spec[i][0][0])
            kin[-1].append(quad_spec[i][0][1])
            for j, item in enumerate(mat):
                mat[j].append(0)
                if quad_spec[i][0][0] == j:
                    mat[j][-1] -= 1
                if quad_spec[i][0][1] == j:
                    mat[j][-1] -= 1
                if quad_spec[i][1][0] == j:
                    mat[j][-1] += 1
                if quad_spec[i][1][1] == j:
                    mat[j][-1] += 1

        for each in quad_species:
            bi_bi_r.append(each[0][0] + ' + ' + each[0][1] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')

    ansatz_str = ''

    rxn_index = 0
    for each in syn_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in deg_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in deg2_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in uni_uni_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in bi_uni_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in uni_bi_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in bi_bi_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    model_str = f"""
functions {{
  vector sys(real t,
             vector y,
             vector theta) {{
      vector[{str(n_species)}] dydt;
      vector[{str(function_count)}] v;
      matrix[{str(n_species)}, {str(function_count)}] S = [
"""
    for each in mat[:-1]:
        model_str += '        [' + str(each[0])
        for item in each[1:]:
            model_str += ', ' + str(item)
    model_str += '],\n'
    model_str += '        [' + str(mat[-1][0])
    for item in mat[-1][1:]:
        model_str += ', ' + str(item)
    model_str += ']\n      ];\n'

    for i in range(len(mat[0])):
        model_str += '      v[' + str(i+1) + '] = theta[' + str(i+1) + ']'
        for each in kin[i]:
            model_str += ' * ' + 'y[' + str(each+1) + ']'
        model_str += ';\n'
    model_str += '      dydt = S * v;\n'
    model_str += '      return dydt;\n'
    model_str += '  }\n'
    model_str += '}'

    model_str += '''
data {
    int N; // Number of observations
    int M; // Number of species
    int obs_idx[M]; // Indices of observed speces
    int D; // Number of possible reactions
    vector[M] y0; // observation at time=0 (used for ode solver)
    real y[N, M]; // Actual data (observations species)
    real ts[N + 1]; // Timepoints
    // horseshoe parameters
    real m0;
    real slab_df;
    real slab_scale;
    // real<lower = 0> tau0;
    // noise model parameters
    // real<lower = 0> noise_sigma;
}
transformed data {
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5 * slab_df;
}
parameters {
    vector<lower = 0>[D] rates_tilde;
    vector<lower = 0>[D] lambda;
    real<lower=0> c2_tilde;
    real<lower=0> tau_tilde;
    real<lower=0> sigma;
}
transformed parameters {
    vector[D] rates;
    vector[D] lambda_tilde;
    vector[M] y_hat[N];
    {
      real tau0 = (m0 / (D - m0)) * (sigma / sqrt(1.0 * N));
      real tau = tau0 * tau_tilde; 
      real c2 = slab_scale2 * c2_tilde;
      lambda_tilde = sqrt( (c2 * square(lambda)) ./ (c2 + square(tau) * square(lambda)) );
      rates = tau * lambda_tilde .* rates_tilde;
    }
     // same as Z(t_j), numerical integration of f(Z(t)) from 0 to t_j
     y_hat = ode_rk45(sys,
                      y0,
                      ts[1],
                      ts[2:],
                      rates);
}
model {
    // horseshoe priors
    rates_tilde ~ normal(0, 1);
    lambda ~ cauchy(0, 1);
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
    tau_tilde ~ cauchy(0, 1);
    // model likelihood
    sigma ~ normal(0, 2);
    for(j in 1:M) {
      y[ ,j] ~ lognormal(log(y_hat[ ,obs_idx[j]]), sigma);
    }

}
generated quantities {
  real y_rep[N, M];
  for(i in 1:N) {
    for(j in 1:M) {
      y_rep[i,j] = lognormal_rng(log(y_hat[i,j]), sigma);
    }
  } 
}
    '''

    return ansatz_str, model_str, function_count
