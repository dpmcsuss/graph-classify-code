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
alg.save=1;


%% incoherent stuff


naive_bayes = find(sum(As,3)>0);
num_inc_edges=1:length(naive_bayes);


% get fisher significant
d_pval = run_get_fisher_pvals(As,constants);
[foo, fis_ind]=sort(d_pval(:),'ascend');

max_ind=find(foo<1-1e-4, 1, 'last' );


% set-up subspaces
subspace.no_sig=fis_ind;
subspace_names  = fieldnames(subspace);
n_subspaces=length(subspace_names);

for m=1:max_ind-1
    disp(m)

    subspace.no_sig=fis_ind(m:max_ind);

    [Lhat_tmp Lsem_tmp] = classify_holdout(As,constants,alg,subspace);
    
    for k=1:n_subspaces
        Lhats.(subspace_names{k})(m)=mean(Lhat_tmp.(subspace_names{k}));
        Lsems.(subspace_names{k})(m)=std(Lsem_tmp.(subspace_names{k}))/sqrt(alg.num_repeats);
    end
    
end

save('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_Lhat_vs_nosig.mat')

%%
% clear, load('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_Lhat_vs_nosig.mat')

prep_plots;

for i=1:n_subspaces
        errorbar(1:max_ind-1,Lhats.(fn{i}),Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{mod(i,3)+1},'linewidth',2);
end

axis('tight')
xlabel('size of removed signal subgraph','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)

alg.save=1;
alg.ext='_Lhat_vs_nosig';
save_plots;
