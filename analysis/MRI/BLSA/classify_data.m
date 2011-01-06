
%% load data


clear, clc, clf
alg.fname='BLSA50';
datadir = '~/Research/data/MRI/BLSA/';
alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';

switch alg.fname
    case 'BLSA79'
        load([alg.datadir '2010_11_07/' alg.fname '_data'])
        alg.s0_trn=[5 10 15 20 30 30];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[5 10 15 20 30 40];     % list of #'s of hold-out samples for class 1
    case 'BLSA42'
        load(['~/Research/data/MRI/BLSA/will/' alg.fname])
        S=length(labels);
        targs=nan(1,S);
        for i=1:S targs(i)=labels(i).gender; end
        alg.s0_trn=[5 10 15 20 20];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[5 10 15 20 25];     % list of #'s of hold-out samples for class 1
        alg.save=0;
    case 'BLSA42B'
        load(['~/Research/data/MRI/BLSA/will/BLSA42_data2'])
        S=length(labels);
        targs=nan(1,S);
        for i=1:S targs(i)=labels(i).gender; end
        alg.s0_trn=[5 10 15 20 20];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[5 10 15 20 25];     % list of #'s of hold-out samples for class 1
        alg.save=0;
    case 'kidney_egg'   % this simulation is easy and should get quite small Lhat
        load('~/Research/data/sims/kidney_egg')
        As=adjacency_matrices;
        targs=class_labels;
        alg.s0_trn=[100 200 500];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[100 200 500];     % list of #'s of hold-out samples for class 1
        alg.save=0;
    case 'kidney_egg2'   % this simulation is easy and should get quite small Lhat
        load('~/Research/data/sims/labeled/kidney_egg')
        As=adjacency_matrices;
        targs=class_labels;
        alg.s0_trn=[4 50 100];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[4 50 100];     % list of #'s of hold-out samples for class 1
        alg.save=0;
    case 'homo_kidney_egg'   % this simulation is easy and should get quite small Lhat
        load(['~/Research/data/sims/unlabeled/' alg.fname])
        As=adjacency_matrices;
        targs=class_labels;
        alg.save=0;
        alg.s0_trn=[10 20 50 100 200];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[10 20 50 100 200];     % list of #'s of hold-out samples for class 1
    case 'hetero'   % this simulation is easy and should get quite small Lhat
        load(['~/Research/data/sims/unlabeled/' alg.fname])
        As=adjacency_matrices;
        targs=class_labels;
        alg.save=0;
    case 'FWQAP'   % this simulation is easy and should get quite small Lhat
        load(['~/Research/data/sims/unlabeled/' alg.fname])
        As=adjacency_matrices;
        targs=class_labels;
        alg.save=0;
    case 'random' % this simulation is impossible and should get Lhat=0.5
        s=80;
        As=round(rand(70,70,s));
        for i=1:s
            As(:,:,i)=triu(As(:,:,i),1);
        end
        targs=[zeros(s/2,1); ones(s/2,1)]';
    case 'BLSA50' % this simulation is impossible and should get Lhat=0.5
        load([alg.postdir alg.fname '_data'])
        s=length(uniq_keeps);
        As=zeros(70,70,s);
        for i=1:s
            AA = uniq_keeps(i).Adj;
            AA(isnan(AA))=0;
            AA(AA<0.2)=0;
            AA(AA>0.2)=1;
            As(:,:,i)=AA;
            targs(i) = uniq_keeps(i).sex;
        end
        alg.s0_trn=[5 10 19 19];     % list of #'s of hold-out samples for class 0
        alg.s1_trn=[5 10 20 25];     % list of #'s of hold-out samples for class 1
end

alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';
alg.save=1;

%% prepare for classification

if length(alg.s0_trn)~=length(alg.s1_trn), error('length of lists of training set sizes must be the same'), end

alg.num_repeats = 100;                   % # of times to repeat using the same # of hold-outs
alg.num_splits  = length(alg.s1_trn);       % # of different hold-outs
alg.num_train_samples=alg.s0_trn+alg.s1_trn;    % total # of training samples

alg.naive_bayes = true;
alg.invars      = false;
alg.ind_edge    = false;

constants = get_constants(As,targs);     % get constants to ease classification code

%% naive bayes prep
if alg.naive_bayes
    subspace.naive_bayes = find(sum(As,3)); % use all edges for which at least 1 edge is present in the data
    subspace.incoherent = [];               % incoherent edge classifier
    subspace.cor_incoherent = [];           % corrected incoherent edge classifier
    subspace.coherent = [];                 % coherent edge classifier
    subspace_names  = fieldnames(subspace);             % pre-allocate space for Lhats
    n_subspaces=length(subspace_names);
    for k=1:n_subspaces
        Lhats.(subspace_names{k})=nan(alg.num_repeats,alg.num_splits);
    end
    
    % test ind-edge using in-sample training data
    [Lhat_in params_in subspace_in] = run_naive_bayes(As,constants,As,constants,subspace);
    disp(Lhat_in)
end

%% graph invariant prep
if alg.invars
    if ~exist('graph_invars','var')
        AAs=As; for i=1:constants.s, AAs(:,:,i)=AAs(:,:,i)'+AAs(:,:,i); end
        graph_invars = get_graph_invariants(AAs);
        save([alg.datadir alg.fname '_invars'],'graph_invars')
    end
    invars.LDA=[];
    invars.QDA=[];
    invars.dLDA=[];
    invars.dQDA=[];
    
    discrimant_names  = fieldnames(invars);             % pre-allocate space for Lhats
    n_discriminants = length(discrimant_names);
    for k=1:n_discriminants
        Lhats.(discrimant_names{k})=nan(alg.num_repeats,alg.num_splits);
    end
else
    graph_invars=[];
end


%% test hold out
Lhats = classify_holdout(As,constants,alg,subspace,graph_invars);


%% make plots
% plot parameter estimates
est_params  = get_params_mle(As,constants);         % estimate parameters from data
plot_model2(est_params,alg);

% plot subspaces
if alg.naive_bayes
    [d_ind deg_ind] = plot_subspaces(est_params.d_pos,subspace,constants.n,alg);
end

% plot Lhats
[Lhat_avg Lhat_sem] = plot_Lhats2(Lhats,alg);                                           % plot misclassification rates

if alg.save
    save([alg.postdir alg.fname '_results'],'est_params','Lhat_avg','Lhat_sem','d_ind','deg_ind','-append')
end