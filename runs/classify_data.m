
clear, clc, clf
%% load data

switch 'BLSA'
    case 'BLSA'
        load('~/Research/data/MRI/BLSA/BLSA79/BLSA79_data')
    case 'kidney_egg'
        load('~/Research/data/sims/kidney_egg')
        As=adjacency_matrices;
        targs=class_labels;
        alg.save=0;
    case 'random'
        s=80;
        As=round(rand(70,70,s));
        for i=1:s
            As(:,:,i)=triu(As(:,:,i),1);
        end
        targs=[zeros(s/2,1); ones(s/2,1)]';
end

alg.save=1;

%% prepare for classification

trn_s0=[5 10 20 30 30];     % list of #'s of hold-out samples for class 0
trn_s1=[5 10 20 30 40];     % list of #'s of hold-out samples for class 1
if length(trn_s0)~=length(trn_s1), error('length of lists of training set sizes must be the same'), end

alg.num_repeats = 10;                   % # of times to repeat using the same # of hold-outs
alg.num_splits  = length(trn_s1);       % # of different hold-outs
alg.num_train_samples=trn_s0+trn_s1;    % total # of training samples

subspace.naive_bayes = find(sum(As,3)); % use all edges for which at least 1 edge is present in the data
subspace.incoherent = [];               % incoherent edge classifier
subspace.cor_incoherent = [];           % corrected incoherent edge classifier
subspace.coherent = [];                 % coherent edge classifier

fn  = fieldnames(subspace);             % pre-allocate space for Lhats
n_subspaces=length(fn);
for k=1:n_subspaces
    Lhats.(fn{k})=nan(alg.num_repeats,alg.num_splits);
end

%% test using in-sample training data
constants = get_constants(As,targs);     % get constants to ease classification code
Lhat_in = run_naive_bayes(As,constants,As,constants,subspace);
disp(Lhat_in)

%% test hold out
for i=1:alg.num_repeats
    disp(i)
    for j=1:length(trn_s0);
        [Atrn Gtrn Atst Gtst] = crossvalprep(As,constants,trn_s0(j),trn_s1(j)); % seperate data into training and testing sets
        
        % ind_edge classifiers
        Lhat_tmp = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace);   % run independent edge classifiers
        for k=1:n_subspaces                                         % store Lhats in a nice way
            Lhats.(fn{k})(i,j)=Lhat_tmp.(fn{k});
        end
        
        
    end
end

if alg.save
    save([alg.datadir 'BLSA79_results'],'Lhats')
end

%% make plots
% plot parameter estimates
est_params  = get_params(As,constants);         % estimate parameters from data
plot_model2(est_params,alg);

% plot subspaces
[d_ind deg_ind] = plot_subspaces(est_params.d_pos,subspace,constants.n,alg);

% plot Lhats
[mean_hat std_hat] = plot_Lhats2(Lhats,alg);                                           % plot misclassification rates

if alg.save
    save([alg.datadir 'BLSA79_results'],'est_params','mean_hat','std_hat','d_ind','deg_ind','-append')
end