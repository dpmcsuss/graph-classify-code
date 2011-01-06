clear, clc, clf
alg.fname='BLSA50';
datadir = '~/Research/data/MRI/BLSA/';
alg.datadir = datadir;
alg.postdir = [alg.datadir 'results/'];
alg.figdir = '~/Research/figs/MRI/BLSA/';

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

constants = get_constants(As,targs);     % get constants to ease classification code


alg.s0_trn = constants.s0-5;     % list of #'s of hold-out samples for class 0
alg.s1_trn = constants.s1-5;     % list of #'s of hold-out samples for class 1
alg.num_repeats = 10000;        % # of times to repeat using the same # of hold-outs

%%

m_max=length(find(sum(As,3)>0));


subspace.fis=[];
subspace.inc=[];
subspace.cor=[];
subspace_names  = fieldnames(subspace);
n_subspaces=length(subspace_names);

subspace3.coh=[];
subspace3.deg=[];
subspace3_names  = fieldnames(subspace3);
n_subspace3=length(subspace_names);

Lhat1=zeros(alg.num_repeats,m_max/5);
Lsem1=zeros(alg.num_repeats,m_max/5);

Lhat2=zeros(alg.num_repeats,1);
Lsem2=zeros(alg.num_repeats,1);

Lhat3=zeros(alg.num_repeats,70);
Lsem3=zeros(alg.num_repeats,70);

parfor i=1:alg.num_repeats
    disp(i)
    [Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,alg.s0_trn,alg.s1_trn); % seperate data into training and testing sets
    
    for m=1:m_max/5
        
        % get parameters
        params = get_params_mle(Atrn,Gtrn);
        d_pval = run_get_fisher_pvals(Atrn,Gtrn);
        
        % fisher significant
        [~, fis_ind]=sort(d_pval(:),'ascend');
        subspace.fis=fis_ind(1:m*5);
        
        % inc
        [~, inc_ind]=sort(params.d_pos(:),'descend');
        subspace.inc=inc_ind(1:m*5);
        
        % z-inc
        [~, cor_ind]=sort(params.d_opt(:),'descend');
        subspace.cor=cor_ind(1:m*5);
        
        [Lhat1(i,m) , ~, ~, Lsem1(i,m)] = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace);   % run independent edge classifiers

    end
    
    subspace2.nb=find(params.d_pos>0);
    [Lhat2(i) , ~, ~, Lsem2(i)] = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace2);
    
    
    for m=1:70
        
        deg=sum(d_pval,1)+sum(d_pval,2)';
        [foo deg_ind] = sort(deg,'ascend');

        ind_max_mat=zeros(70);
        ind_max_mat(deg_ind(1:m),deg_ind(1:m))=1;
        ind_max_mat=triu(ind_max_mat,+1);
        subspace3.deg=find(ind_max_mat);
        subspace3.coherent = m;

        [Lhat3(i,m) , ~, ~, Lsem3(i,m)] = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace3);   % run independent edge classifiers
        
    end
end



%%
fn=subspace_names;
fs=12;

prep_plots;

for m=1:70, blowup(m)=choose(m,2); end

figure(1), clf, hold on
for i=1:n_subspaces
    if ~strcmp(fn{i},'coherent') && ~strcmp(fn{i},'deg')
        plot(1:m_max,Lhats.(fn{i}),... Lsems.(fn{i}),...
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

alg.save=0;
alg.ext='_Lhat_vs_m_in';
save_plots;