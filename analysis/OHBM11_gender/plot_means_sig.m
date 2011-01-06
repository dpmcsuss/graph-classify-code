clear, clc


alg.save=1;
alg.figdir = '/Users/jovo/Research/figs/MRI/BLSA/';
alg.fname = 'BLSA50';


%% plot model and significance
clc
figure(1), clf,
fs=10;

load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')


alg.save=1;
alg.figdir = '/Users/jovo/Research/figs/MRI/BLSA/';
alg.fname = 'BLSA50';


params = get_params_mle(As,constants);

emax=max(max([params.E0(:) params.E1(:) abs(params.E0(:)-params.E1(:))]));

subplot(131)
image(64*params.E0/emax)
colormap('gray')
title('female','fontsize',fs)
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

subplot(132)
image(64*params.E1/emax)
title('male','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

subplot(133)
imagesc(1-d_pval)
title('significance of difference','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])



if alg.save
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_params'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end

%% plot Lhat vs signal subgraph size

clear, clc
load('/Users/jovo/Research/data/MRI/BLSA/results/BLSA50_Lhat_200_subspaces.mat')
alg.save=1;
plot_Lhat_vs_signal_subgraph_size

%% make table of m most significant edges

clear, clc,
load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')

[foo IX] = sort(d_pval(:));

[I J] = ind2sub([70 70],IX);

CorticalLabels = importdata('/Users/jovo/Research/data/MRI/BLSA/misc/IACL-corticallabeldescriptions-1110.txt',' ');


clc

clear orderedregions
m=find(foo<0.999999999999, 1, 'last' );
for i=1:m
    
    if I(i)<36, ii=I(i)+2; else ii=I(i)+3; end
    if J(i)<36, jj=J(i)+2; else jj=J(i)+3; end
    
    if ii<11, ii_start=7; elseif ii<39, ii_start=8; else ii_start=9; end
    if jj<11, jj_start=7; elseif jj<39, jj_start=8; else jj_start=9; end
    
    maxlen=19;
    label1=CorticalLabels{ii}(ii_start:end-1);
    if length(label1)>maxlen, label1=label1(1:maxlen); end
    label2=CorticalLabels{jj}(jj_start:end-1);
    if length(label2)>maxlen, label2=label2(1:maxlen); end
    
    orderedregions{i,1}=foo(i);
    orderedregions{i,2}=label1;
    orderedregions{i,3}=label2;
    
end

%% remove good edges

clear, clc,
load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')

alg.save=0;
alg.num_repeats=100;
alg.s0_trn=15;
alg.s1_trn=15;
[pvals ind]=sort(d_pval(:),'ascend');

k=0;
for m=1:500
    k=k+1;
    subspace.manual=ind(m:end);
    [Lhat Lsem] = classify_holdout(As,constants,alg,subspace);
    
    Lhats.manual(k)=mean(Lhat.manual);
    Lsems.manual(k)=std(Lsem.manual)/sqrt(alg.num_repeats);
end


% alg.num_repeats=1000;
% subspace.manual=ind(50:end)
% [LhatX LsemX] = classify_holdout(As,constants,alg,subspace);
% mean(LhatX.manual)


alg.save=0;
alg.figdir = '/Users/jovo/Research/figs/MRI/BLSA/';
% alg.postdir = '/Users/jovo/Research/data/MRI/BLSA/results/';
alg.fname = 'BLSA50';
if alg.save
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_removed_subspace'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end


%% Lhat vs. m

clear, clc
load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')


alg.save=1;
alg.num_repeats=100;
alg.s0_trn=18;
alg.s1_trn=18;
alg.postdir='/Users/jovo/Research/data/MRI/BLSA/results/';
alg.fname='Lhat_v_m';

naive_bayes=find(sum(As,3)>0);


[foo fis_ind]=sort(d_pval(:),'ascend');

params = get_params_mle(As,constants);

[foo inc_ind]=sort(params.d_pos(:),'ascend');
[foo cor_ind]=sort(params.d_opt(:),'ascend');


subspace.fis=[];
subspace.inc=[];
subspace.cor=[];

subspace_names  = fieldnames(subspace);
n_subspaces=length(subspace_names);

for m=1:length(naive_bayes)
    
    disp(m)
    subspace.fis=fis_ind(1:m);
    subspace.inc=inc_ind(1:m);
    subspace.cor=cor_ind(1:m);
    [Lhat Lsem] = classify_holdout(As,constants,alg,subspace);
    
    for k=1:n_subspaces
        Lhats.(subspace_names{k})(m)=mean(Lhat.(subspace_names{k}));
        Lsems.(subspace_names{k})(m)=std(Lsem.(subspace_names{k}))/sqrt(alg.num_repeats);
    end
    
end


%%
deg=sum(d_pval)+sum(d_pval');
[foo deg_ind] = sort(deg,'ascend');

clear subspace Lhat Lsem
for m=5:length(deg_ind)
    
    disp(m)
    
    ind_max_mat=zeros(70);
    ind_max_mat(deg_ind(1:m),deg_ind(1:m))=1;
    ind_max_mat=triu(ind_max_mat,+1);
    subspace.coh=find(ind_max_mat);
    
    [Lhat Lsem] = classify_holdout(As,constants,alg,subspace);
    
    Lhats.coh(m)=mean(Lhat.coh);
    Lsems.coh(m)=std(Lsem.coh)/sqrt(alg.num_repeats);
    
end











