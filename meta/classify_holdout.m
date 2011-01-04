function [Lhats Lsems alg] = classify_holdout(As,constants,alg,subspace,graph_invars)
% this function computes trains and classifies using held-out data
% unbalanced training data is permitted
% number test samples can be a function of number of training samples
%
% INPUT:
%   As:         3rd-order tensor of adjacency_matrices
%   constants:  structure containing labels as well as other info
%   alg:        specifies which algorithms to use, and which parameters for each alg
%   subspace:   (optional) specifies which subspaces to use if naive_bayes
%               classifier is used
%   graph_invars: (optional) specifies the graph invariants to be used
%
% OUTPUT:
%   Lhat: estimate misclassification rates

if isempty(subspace),           nb_do=0;                else nb_do=1; end
if nargin==5,                   invars_do=1;            else invars_do=0; end
if isfield(alg,'save'),         save_do=alg.save;       else save_do=0; end
if ~isfield(alg,'num_repeats'), alg.num_repeats=10;     end
if ~isfield(alg,'s0_trn'),      alg.s0_trn=round(constants.s0*0.9); end
if ~isfield(alg,'s1_trn'),      alg.s1_trn=round(constants.s1*0.9); end
if ~isfield(alg,'trn_s'),       alg.trn_s=alg.s0_trn+alg.s1_trn; end


if nb_do
    subspace_names  = fieldnames(subspace);             % pre-allocate space for Lhats
    n_subspaces=length(subspace_names);
end

for i=1:alg.num_repeats
    if mod(i,100)==0, disp(['repeat # ' num2str(i)]), end
    for j=1:length(alg.s0_trn);
        [Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,alg.s0_trn(j),alg.s1_trn(j)); % seperate data into training and testing sets
        
        % naive_bayes classifiers
        if nb_do
            [Lhat_tmp foo foo Lsem_tmp] = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace);   % run independent edge classifiers
            for k=1:n_subspaces                                         % store Lhats in a nice way
                Lhats.(subspace_names{k})(i,j)=Lhat_tmp.(subspace_names{k});
                Lsems.(subspace_names{k})(i,j)=Lsem_tmp.(subspace_names{k});
            end
        end
        
        % graph invariants classifers
        if invars_do
            QDA_params = get_discriminant_params(graph_invars,inds);
            Lhat_tmp = discriminant_classifiers(graph_invars,targs,QDA_params,invars);
            for k=1:n_discriminants                                         % store Lhats in a nice way
                Lhats.(discrimant_names{k})(i,j)=Lhat_tmp.(discrimant_names{k});
            end
        end
    end
end

Lhats=orderfields(Lhats);
Lsems=orderfields(Lsems);

if save_do
    save([alg.postdir alg.fname '_results'],'Lhats')
end
