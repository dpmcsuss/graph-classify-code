% %% Pipeline as outlined below by Jovo (via e-mail) 
%
% 1) for each graph, compute the scaled sign positive graph laplacian of
% each graph (call it L_i).  one obtains this by filling in the diagonal of
% each adjacency matrix by the degree of that vertex, divided by n-1 (this
% is the expected weight of the self loop).
% 
% 2) compute an svd of L_i for all i.  look at the scree plots.
% presumably, a relatively small number of eigenvalues should suffice to
% capture most of the variance. 
% 
% 3) build a classifier using only the first d singular vectors.  we
% haven't really thought much about this, but svd implicitly assumes
% gaussian noise, so, the optimal classifier under this assumption is a
% discriminant classifier.  the code i have in the repo for discriminant
% classifiers is still in progress, but if you want to fix it up, that'd be
% fine.  report results for the discriminant classifiers (LDA and QDA, but
% with diagonal covariance matrices and whitened covariance matrices).
% 
% 4) cluster the singular vectors.  do this by initializing with a
% hierarchical clustering algorithm, and then use that to seed a k-means
% algorithm.  you can do this for a variety of values of k and d.  choose
% the one with the best BIC. marchette will send details for computing BIC
% values here, and for implementing the "k-means" algorithm, which is, in
% fact, more like a constrained finite gaussian mixture model algorithm.
% 
% 5) given a particular choice of k and d, build a classifier.  here, each
% class is a mixture of gaussians, so a naive discriminant classifier is
% insufficient.  instead, one can implement a bayes plugin classifier. so,
% given the fit of the mixture distributions for the two classes, compute
% whether a test graph is closer to class 0 or class 1 distribution.

%% Setup - Change to fit your desires

% this is just where i load stuff from my computer
load /users/dsussman/documents/MATLAB/Brain_graphs/BLSA79_data.mat;

% you can set up stuff however you want here

%% Step 1 - Compute Laplacians
Ls = graph_laplacian(As);

%% Step 2 - Compute SVDs and make scree plots
[Us, Ss, Vs] = graph_svd(Ls);
sz =  size(Ss);
singVals =  zeros(sz(1),sz(3)); % col vectors of singular vals
for k=1:sz(3)
    singVals(:,k) = diag(squeeze(Ss(:,:,k)));
end


% show scree plots of for all the graphs, red is targs is true blue is
% targs is false
figure(401);
hold all;
plot( singVals(:,targs==1), 'r');
plot( singVals(:,targs==0), 'b');

%% Step 3 - Classify using first d singular values
% needed to add the repo to the path
d=5;
inds = struct('ytrn', 1:length(targs),...
              'y0trn', find(targs == 0),...
              'y1trn', find(targs == 1));
discrim = struct('dLDA',1);
params = get_discriminant_params(singVals(1:d,:),inds,discrim);

[Lhat Lsem] = discriminant_classifiers(singVals(1:d,:),targs,params,discrim);