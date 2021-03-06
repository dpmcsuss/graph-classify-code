function P = get_subspace_params(adjacency_matrices,constants,subspace)
% this funtion gets the parameters necessary for independent edge algorithms

% NOTE THAT THIS IS ENSURES WE DON'T GET ZEROS
eps = 1/(10*constants.s);                        % to deal with 0's and 1's

% estimated class 0 edge probabilities
P.E0    = mean(adjacency_matrices(:,:,constants.y0),3);
P.E0(P.E0==0) = eps;
P.E0(P.E0==1) = 1-eps;

% estimated class 0 edge probabilities
P.E1    = mean(adjacency_matrices(:,:,constants.y1),3);
P.E1(P.E1==0) = eps;
P.E1(P.E1==1) = 1-eps;

% pre-compute constants for bernoulli distribution
P.lnE0  = log(P.E0);
P.ln1E0 = log(1-P.E0);
P.lnE1  = log(P.E1);
P.ln1E1 = log(1-P.E1);

% log-priors
P.lnprior0 = log(constants.s0/constants.s);
P.lnprior1 = log(constants.s1/constants.s);