function P = get_params_map(adjacency_matrices,constants)
% this funtion gets the parameters necessary for independent edge algorithms

% set prior to be flat
alpha = 1;
beta = 1;


% get posteriors
sum0 = sum(adjacency_matrices(:,:,constants.y0)); 
sum1 = sum(adjacency_matrices(:,:,constants.y1)); 

P.alpha0    = alpha + sum0;
P.alpha1    = alpha + sum1;

P.beta0     = beta + constants.s0 - sum0;
P.beta1     = beta + constants.s1 - sum1;

P.mode0     = (P.alpha0-1)/(P.alpha0+P.beta0-2);
P.mode1     = (P.alpha1-1)/(P.alpha1+P.beta1-2);


% various measures to compute difference matrix
P.d_KLbeta = KLbeta(P);

% pre-compute constants for bernoulli distribution
P.lnE0  = log(P.E0);
P.ln1E0 = log(1-P.E0);
P.lnE1  = log(P.E1);
P.ln1E1 = log(1-P.E1);

% log-priors
P.lnprior0 = log(constants.s0/constants.s);
P.lnprior1 = log(constants.s1/constants.s);


end

function KL = KLbeta(P)

a0=P.alpha0;
a1=P.alpha1;
b0=P.beta0;
b1=P.beta1;

KL = log(beta(a1,b1)/beta(a0,b0))-(a1-a0)*digraph(a0)-(b1-b0)*digraph(b0)+(a1-a0+b1-b0)*digraph(a0+b0);


end