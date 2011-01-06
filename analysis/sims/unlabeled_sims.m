%% simulate independent edge models and classify using or not using the vertex names
clear; clc

alg.datadir = '~/Research/data/graph_sims/unlabeled/';
alg.figdir  = '~/Research/figs/graph_sims/unlabeled/';
alg.fname   = 'homo_kidney_egg';    % different names will generate different simulations

n   = 10;                           % # of vertices
S   = 500;                          % # of samples

alg.save    = 0;                    % whether to save/print results
alg.names   = [{'LAP'}; {'QAP'}];   % which algorithms to run
alg.truth_start = 1;                % start LAP and QAP at truth

alg.QAP_max_iters   = 2;            % max # of iterations when using QAP approx
alg.QAP_init        = eye(n);       % starting value for QAP (see sfw comments for explanation)

switch alg.fname
    case 'homo_kidney_egg'
        
        p   = 0.5;      % prob of connection for kidney
        q0  = 0.25;     % prob of connection for egg
        q1  = 0.75;     % prob of connection for egg
        egg = 1:3;      % vertices in egg
        
        E0=p*ones(n);   % params in class 0
        E0(egg,egg)=q0; % egg params in class 0
        
        E1=p*ones(n);   % params in class 1
        E1(egg,egg)=q1; % egg params in class 1
                
        params.n=n; params.p=p; params.q0=q0; params.q1=q1; params.egg=egg; params.S=S; params.E0=E0; params.E1=E1;
        
    case 'hetero'
         
        E0=rand(n);     % params in class 0
        E1=rand(n);     % params in class 1
                
        params.n=n; params.S=S; params.E0=E0; params.E1=E1;
        
    case 'hard_hetero'
         
        E0=rand(n);         % params in class 0
        E1=E0+randn(n)*.1;  % params in class 1
        E1(E1>=1)=1-1e-3;   % don't let prob be >1
        E1(E1<=0)=1e-3;     % or <0
                
        params.n=n; params.S=S; params.E0=E0; params.E1=E1;
end

S0 = S/2; % # of samples in class 0
S1 = S/2; % # of samples in class 1

A0 = repmat(E0,[1 1 S0]) > rand(n,n,S0);    % class 0 samples
A1 = repmat(E1,[1 1 S1]) > rand(n,n,S1);    % class 1 samples

adjacency_matrices(:,:,1:S0)=A0;            % create adjacency_matrices to store all sampled graphs
adjacency_matrices(:,:,S0+1:S)=A1;          % add class 1 samples
class_labels=[zeros(1,S0) ones(1,S1)];      % vector of class labels

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')

%% performance tests

trn_ind = [1:S0/2,      S0+1:S0+S1/2];  % let half of the samples from each class be training
tst_ind = [S0/2+1:S0,   S0+S1/2+1:S];   % the other half is for testing

ytst = class_labels(tst_ind);           % class labes for test data

% parameters for naive bayes classifiers
P.lnE0  = log(params.E0);
P.ln1E0 = log(1-params.E0);
P.lnE1  = log(params.E1);
P.ln1E1 = log(1-params.E1);

P.lnprior0 = log(S0/S);
P.lnprior1 = log(S1/S);

%% when not trying to solve GIP

% make data unlabeled
adj_mat=0*adjacency_matrices;
for i=S
    q=randperm(n);
    A=adjacency_matrices(:,:,i);
    adj_mat(:,:,i)=A(q,q);
end

Atst=adj_mat(:,:,tst_ind);

[Lhat Lvar] = naive_bayes_classify(Atst,ytst,P);
Lhats.rand  = Lhat.all;
Lvars.rand  = Lvar.all;

%% performance using true parameters and labels

Atst=adjacency_matrices(:,:,tst_ind);

[Lhat Lvar] = naive_bayes_classify(Atst,ytst,P);
Lhats.star  = Lhat.all;
Lvars.star  = Lvar.all;

%% test using hold-out training data

Atrn=adjacency_matrices(:,:,trn_ind);

[LAP QAP] = classify_unlabeled_graphs(Atrn,Atst,ytst,P,alg);

save([alg.datadir alg.fname '_results'],'adjacency_matrices','class_labels','params','alg','Lhat','LAP','QAP')

%% plot model

figure(2), clf
fs=10;  % fontsize

emax=max(max([params.E0(:) params.E1(:) abs(params.E0(:)-params.E1(:))]));

% class 0
subplot(131)
image(60*params.E0/emax)
colormap('gray')
title('class 0 mean','fontsize',fs)
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% class 1
subplot(132)
image(60*params.E1/emax)
title('class 1 mean','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

% difference
subplot(133)
image(60*abs(params.E0-params.E1)/emax)
colormap('gray')
title('difference','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
colorbar('Position',[.925 .3 .02 .4],'YTick',[0 15 30 45 60],'YTickLabel',[0 25 50 75 100],'fontsize',fs)


if alg.save
    wh=[6 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_model'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end

%% plot Lhat & errorbars

figure(3), clf, hold all
ms=16;

% rand
errorbar(0,Lhats.rand,Lvars.rand,'g','linewidth',2,'Marker','.','Markersize',ms)

% LAP
if strcmp(alg.names(1),'LAP'),
    errorbar(1.1,LAP.Lhat,LAP.Lvar,'k','linewidth',2,'Marker','.','Markersize',ms)
end

% QAP
if strcmp(alg.names(2),'QAP'),
    errorbar(1:alg.QAP_max_iters+.2,QAP.Lhat,QAP.Lvar,'linewidth',2,'Marker','.','Markersize',ms)
end

% L*
errorbar(alg.QAP_max_iters+1,Lhats.star,Lvars.star,'r','linewidth',2,'Marker','.','Markersize',ms)

legend('rand','LAP','QAP','L^*')
axis([-0.5 alg.QAP_max_iters+1.5 0 0.5])

ylabel('$\hat{L}$','interp','latex','Rotation',0)
xlabel('# of interations','interp','none')

if alg.save
    wh=[6 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_Lhats'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end


%% plot objective functions for QAP

figure(4), clf, hold all

errorbar(0:alg.QAP_max_iters,QAP.obj0_avg,QAP.obj0_var,'k','linewidth',2)
errorbar(0:alg.QAP_max_iters,QAP.obj1_avg,QAP.obj1_var,'r','linewidth',2)

legend('class 0','class 1')

ylabel('objective function')
xlabel('# of interations','interp','none')

if alg.save
    wh=[6 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_obj'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
