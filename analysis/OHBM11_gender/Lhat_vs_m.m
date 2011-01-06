clear, clc, clf
alg.fname='BLSA50';
datadir = '~/Research/data/MRI/BLSA/';
alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';

load([alg.postdir alg.fname '_data'])
s=length(uniq_keeps);
As=nan(70,70,s);
targs=nan(70);
for i=1:s
    AA = uniq_keeps(i).Adj;
    AA(isnan(AA))=0;
    AA(AA<0.2)=0;
    AA(AA>0.2)=1;
    As(:,:,i)=AA;
    targs(i) = uniq_keeps(i).sex;
end

constants = get_constants(As,targs);     % get constants to ease classification code

alg.s0_trn = constants.s0-5;     % list of #'s of hold-out samples for class 0
alg.s1_trn = constants.s1-5;     % list of #'s of hold-out samples for class 1
alg.num_repeats = 1000;        % # of times to repeat using the same # of hold-outs



%% incoherent stuff

clear subspace

naive_bayes = find(sum(As,3)>0);
num_inc_edges=1:length(naive_bayes);


% get fisher significant
d_pval = run_get_fisher_pvals(As,constants);
[~, fis_ind]=sort(d_pval(:),'ascend');

% get pos and opt rankings
params = get_params_mle(As,constants);
[~, inc_ind]=sort(params.d_pos(:),'descend');
[~, cor_ind]=sort(params.d_opt(:),'descend');

% choose subspaces
subspace.inc = [];
subspace.fis = [];
subspace.cor = [];

% for update Lhats
subspace_names  = fieldnames(subspace);
n_subspaces=length(subspace_names);

for m=1:length(num_inc_edges)
    disp(num_inc_edges(m))
    
    subspace.inc=inc_ind(1:m);
    subspace.cor=cor_ind(1:m);
    subspace.fis=fis_ind(1:m);
    
    [Lhat_tmp Lsem_tmp] = classify_holdout(As,constants,alg,subspace);
    
    for k=1:n_subspaces
        Lhats.(subspace_names{k})(m)=mean(Lhat_tmp.(subspace_names{k}));
        Lsems.(subspace_names{k})(m)=std(Lsem_tmp.(subspace_names{k}))/sqrt(alg.num_repeats);
    end
    
end

%%

clear subspace

num_coh_vertices=1:70;

deg=sum(d_pval)+sum(d_pval');
[foo deg_ind] = sort(deg,'ascend');

subspace.coherent = [];
subspace.deg = [];

subspace_names  = fieldnames(subspace);
n_subspaces=length(subspace_names);

for m=1:length(num_coh_vertices)
    disp(num_coh_vertices(m))
    
    subspace.coherent = num_coh_vertices(m);
    
    ind_max_mat=zeros(70);
    ind_max_mat(deg_ind(1:m),deg_ind(1:m))=1;
    ind_max_mat=triu(ind_max_mat,+1);
    subspace.deg=find(ind_max_mat);
    
    [Lhat_tmp Lsem_tmp] = classify_holdout(As,constants,alg,subspace);
    
    for k=1:n_subspaces
        Lhats.(subspace_names{k})(m)=mean(Lhat_tmp.(subspace_names{k}));
        Lsems.(subspace_names{k})(m)=std(Lsem_tmp.(subspace_names{k}))/sqrt(alg.num_repeats);
    end
    
end

%%

clear, load('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_Lhat_vs_m.mat')

prep_plots;

for m=1:70, blowup(m)=choose(m,2); end

figure(1), clf, hold on
for i=1:n_subspaces
    if ~strcmp(fn{i},'coherent') && ~strcmp(fn{i},'deg')
        plot(num_inc_edges,Lhats.(fn{i}),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{i},'linewidth',2);
    else %if strcmp(fn{i},'coherent')
        plot(blowup,Lhats.(fn{i}),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{mod(i,3)+1},'linewidth',2);
    end
end

axis([1 length(naive_bayes) 0 0.51])
legend('inc','fisher','norm inc','coh_d_e_g','coh_f_i_s','location','best')
set(gca,'XScale','log')
xlabel('size of signal subgraph','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)

alg.save=1;
alg.ext='_Lhat_vs_m';
save_plots;