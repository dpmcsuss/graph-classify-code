function [Lhats xval] = classify_par(As,constants,xval,subspace,preds)
% this function computes trains, classifies and estimates performance of
% each of the various classifiers using cross-validated data.
%
% INPUT:
%   As:         array of s adjacency_matrices, each n x n
%   constants:  structure containing labels as well as other info
%   xval:       specifies at a minimum: s0_trn and s1_trn
%   subspace:   (optional) specifies which subspaces to use for plugin classifiers
%   preds:      (optional) specifies predictive variables to be used for discriminant classifiers
%
% OUTPUT:
%   Lhats:      estimated misclassification rates
%   xval:       xval structure after changing it

%% ensure that xval fields are set
% if ~isfield(xval,'s0_trn'),
%     disp(['s0 = ' num2str(constants.s0)])
%     xval.s0_trn = input('\s0_trn [positive integer]: ');
% end
% if ~isfield(xval,'s1_trn'),
%     disp(['s1 = ' num2str(constants.s1)])
%     xval.s1_trn = input('\s1_trn [positive integer]: ');
% end

% if num_iters is not specified, do leave-one-out cross-validation
if ~isfield(xval,'num_iters'),
    xval.num_iters = constants.s; %choose(constants.s0,xval.s0_trn)*choose(constants.s1,xval.s1_trn);
end



%% run classifiers
Lstats_tmp(1:xval.num_iters)=struct('avg',nan,'sem',nan);
for i=1:xval.num_iters
    if mod(i,100)==0, disp(['repeat # ' num2str(i)]), end
    [Atrn Gtrn Atst Gtst inds] = crossvalprep3(As,constants,xval,i); % seperate data into training and testing sets
    
    % naive_bayes classifiers
    if ~isempty(subspace)
        Lstats_tmp(i)= run_plugin_classifier(Atrn,Gtrn,Atst,Gtst,subspace);   % run independent edge classifiers
    end
    
    % graph invariants classifers
%     if nargin==5
%         discrim_params = get_discriminant_params(preds,inds);
%         Lhat.discrim(i) = discriminant_classifiers(preds,targs,discrim_params,invars);
%     end
end

%% store Lstats
Lhats=nan(xval.num_iters,length(subspace));
for i=1:xval.num_iters
    Lhats(i,:)=Lstats_tmp(i).avg;
end