%% simulate kidney-egg problem
clear; clc

alg.datadir = '~/Research/data/graph_sims/unlabeled/';
alg.figdir  = '~/Research/figs/graph_sims/unlabeled/';
alg.save    = 1; % whether to save/print results

sim_name = 'hetero'; %'homo_kidney_egg'
alg.fname   = sim_name;

switch sim_name
    case 'homo_kidney_egg'
        
        n   = 10;   % # of vertices
        p   = 0.5;  % prob of connection for kidney
        q0  = 0.25; % prob of connection for egg
        q1  = 0.75; % prob of connection for egg
        egg = 1:3;  % vertices in egg
        S   = 20000;  % # of samples
        
        E0=p*ones(n);   % params in class 0
        E0(egg,egg)=q0; % egg params in class 0
        
        E1=p*ones(n);   % params in class 1
        E1(egg,egg)=q1; % egg params in class 1
        
        A0 = repmat(E0,[1 1 S/2]) > rand(n,n,S/2); % class 0 samples
        A1 = repmat(E1,[1 1 S/2]) > rand(n,n,S/2); % class 1 samples
        
        siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
        siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
        adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
        adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
        class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels
        
        params.n=n; params.p=p; params.q0=q0; params.q1=q1; params.egg=egg; params.S=S; params.E0=E0; params.E1=E1;
        
    case 'hetero'
        
        n   = 10;   % # of vertices
        S   = 20000;  % # of samples
        
        E0=rand(n);   % params in class 0
        E1=rand(n);   % params in class 1
        
        A0 = repmat(E0,[1 1 S/2]) > rand(n,n,S/2); % class 0 samples
        A1 = repmat(E1,[1 1 S/2]) > rand(n,n,S/2); % class 1 samples
        
        siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
        siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
        adjacency_matrices(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
        adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
        class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels
        
        params.n=n; params.S=S; params.E0=E0; params.E1=E1;
        
end

save([alg.datadir alg.fname],'adjacency_matrices','class_labels','params','alg')

%% knobs to tweak in classifiers

alg.names = [{'LAP'}; {'QAP'}]; % which algorithms to run
alg.QAP_max_iters = 15;         % max # of iterations when using QAP approx
alg.nb_ind = 1:n^2;             % use all edges in labeled classifiers

%% performance tests

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

trn_ind = [1:constants.s0/2,                constants.s0+1:constants.s0+constants.s1/2];
tst_ind = [constants.s0/2+1:constants.s0,   constants.s0+constants.s1/2+1:constants.s];

ytrn=constants.ys(trn_ind);
ytst=constants.ys(tst_ind);

adj_mat=0*adjacency_matrices;
% make data unlabeled
for i=constants.s
    q=randperm(constants.n);
    A=adjacency_matrices(:,:,i);
    adj_mat(:,:,i)=A(q,q);
end

%% when not trying to solve GIP

Atrn=adj_mat(:,:,trn_ind);
Atst=adj_mat(:,:,tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

[Lhat,~,Lvar] = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
Lhats.rand = Lhat.nb;
Lvars.rand = Lvar.nb;

%% performance using true parameters and labels

Atrn=adjacency_matrices(:,:,trn_ind);
Atst=adjacency_matrices(:,:,tst_ind);

Gtrn=get_constants(Atrn,ytrn);
Gtst=get_constants(Atst,ytst);

[Lhat,~,Lvar] = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst,params);
Lhats.star = Lhat.nb;
Lvars.star = Lvar.nb;

%% test using hold-out training data
[LAP QAP] = classify_unlabeled_graphs(adjacency_matrices,class_labels,alg,params);

save([alg.datadir alg.fname '_results'],'adjacency_matrices','class_labels','params','alg','Lhat','LAP','QAP')

%% make plots

plot_model(params,alg)                              % plot params and estimated params

%% plot Lhat & errorbars

figure(3), clf, hold all
ms=16;
errorbar(0,Lhats.rand,Lvars.rand,'g','linewidth',2,'Marker','.','Markersize',ms)

if strcmp(alg.names(1),'LAP'),
    errorbar(1.1,LAP.Lhat,LAP.Lvar,'k','linewidth',2,'Marker','.','Markersize',ms)
end

if strcmp(alg.names(2),'QAP'),
    errorbar([1:alg.QAP_max_iters]+.2,QAP.Lhat,QAP.Lvar,'linewidth',2,'Marker','.','Markersize',ms)
end

errorbar(alg.QAP_max_iters+1,Lhats.star,Lvars.star,'r','linewidth',2,'Marker','.','Markersize',ms)

legend('rand','LAP','QAP','L^*')
axis([-0.5 alg.QAP_max_iters+1.5 0 1])

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

errorbar(QAP.obj0_avg,QAP.obj0_var,'k','linewidth',2)
errorbar(QAP.obj1_avg,QAP.obj1_var,'r','linewidth',2)

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
