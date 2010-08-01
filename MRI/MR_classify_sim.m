clear, clc, clf
% rand('seed',12346);

load('~/Research/necog/data/BLSA/will/BLSA32')

As(As<=0.2)=0;
As(As>0.2)=1;

G = get_constants(As,ys);
P = Bin_get_params(As,G);

n       = G.n;          % # vertices
eps     = 0.05;         % don't use zeros and ones

% make E0
E0          = P.E0; 
% E0(E0<=0)   = eps;
% E0(E0>=1)   = 1-eps;

% make E1
E1          = P.E1; 
% E1(E1<=0)   = eps;
% E1(E1>=1)   = 1-eps;

s_trn   = 200;                                              % # training samples
A0      = repmat(E0,[1 1 s_trn/2]) > rand(n,n,s_trn/2);     % class 0 training samples
A1      = repmat(E1,[1 1 s_trn/2]) > rand(n,n,s_trn/2);     % class 1
A_trn   = zeros(n,n,s_trn);                                 % pre-allocate memory
A_trn(:,:,1:2:s_trn) = A0;                                  % training adjacency matrices
A_trn(:,:,2:2:s_trn) = A1;
y_trn   = zeros(1,s_trn);                                   % training labels
y_trn(2:2:s_trn) = 1;
G_trn   = get_constants(A_trn,y_trn);                       % get constants to ease classification code


s_tst   = 500;                                              % # of test samples
A0      = repmat(E0,[1 1 s_tst/2]) > rand(n,n,s_tst/2);
A1      = repmat(E1,[1 1 s_tst/2]) > rand(n,n,s_tst/2);
A_tst   = zeros(n,n,s_tst);
A_tst(:,:,1:2:s_tst) = A0;
A_tst(:,:,2:2:s_tst) = A1;
y_tst   = zeros(1,s_tst);
y_tst(2:2:s_tst) = 1;
G_tst   = get_constants(A_tst,y_tst);

datadir = '~/Research/necog/data/BLSA/';
figdir  = '~/Research/necog/figs/MR/';
fname   = 'MR32_sim';
save_stuff = 1;

clear A0 A1 
if save_stuff, save([datadir fname]); end


%% in sample

Z.nb    = 1;
Z.inc   = 1; Z.Ninc=G.n;
Z.max   = 1; Z.Nmax=round(G.n/10);
Z.case='in';

[Lhatin ind P yhatin] = graph_classify_ie(A_trn,G_trn,Z);

%% hold-out

clear Z
Z.nb    = 1;
Z.inc   = 1; Z.Ninc = G.n;
Z.max   = 1; Z.Nmax = G.n/10;
Z.case  = 'ho';

Ntrain=[10 25 50 75 100];

Lhats.nb=zeros(length(Ntrain),1);
Lhats.inc=zeros(length(Ntrain),1);
Lhats.max=zeros(length(Ntrain),1);

for i=1:length(Ntrain)
                
        Atemp   = A_trn(:,:,1:Ntrain(i)*2);
        ytemp   = y_trn(1:Ntrain(i)*2); 
        Gtemp   = get_constants(Atemp,ytemp);
        
        [Lhat Lvar ind P yhat] = graph_classify_ie(Atemp,Gtemp,Z,A_tst,G_tst);

        Lhats.nb(i)=Lhat.nb;
        Lhats.inc(i)=Lhat.inc;
        Lhats.max(i)=Lhat.max;

        Lvars.nb(i)=Lvar.nb;
        Lvars.inc(i)=Lvar.inc;
        Lvars.max(i)=Lvar.max;

end


if save_stuff, save([datadir fname]); end

%% plot sim

figure(1), clf
nrows=1;
ncols=3;

emax=max(max([E0(:) E1(:) abs(E0(:)-E1(:))]));
subplot(nrows,ncols,1), cla
image(60*E0/emax)
colormap('gray')
title('E_0')
xlabel('edges')
ylabel('edges')

subplot(nrows,ncols,2), cla
image(60*E1/emax)
title('E_1')

subplot(nrows,ncols,3), cla
image(60*abs(E0-E1)/emax)
colormap('gray')
title('\Delta')

if save_stuff
    wh=[6 1.8];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir fname '_params'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% plot misclassification and edge detection rates
figure(888), clf, hold all

ls1='-'; ls2='--'; ls3=':'; ls4='.';
color1=0*[1 1 1]; color2=0.75*[1 1 1]; color3=0.7*[1 1 1]; color4=0.25*[1 1 1];
lw=2;
fs=12;
xticks=Ntrain;
clear Ncorrect ind_IE_mat ind_max_mat

h(1)=errorbar(Ntrain,Lhats.nb,Lvars.nb,'color',color1,'linestyle',ls1,'linewidth',lw);
h(2)=errorbar(Ntrain+1,Lhats.inc,Lvars.inc,'color',color2,'linestyle',ls2,'linewidth',lw);
h(3)=errorbar(Ntrain+2,Lhats.max,Lvars.max,'color',color3,'linestyle',ls3,'linewidth',lw);
% h(4)=plot(Ntrain,Lhats.tru,'color',color4,'linestyle',ls4,'linewidth',lw);
axis([0 max(Ntrain) 0 0.5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)
legend(h,'nb','inc','max','location','best')
set(gca,'XTick',Ntrain,'XTickLabel',2*Ntrain)



if save_stuff
    wh=[5 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'MR_rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% plot example recovery image

figure(11), clf
nrows=3;
ncols=1;
i=0;

load('~/Research/necog/data/BLSA/MR_14')
G = get_constants(As,ys);
P = Bin_get_params(As,G);

% true signal
i=i+1; subplot(nrows,ncols,i), cla
imagesc(P.del)
colormap('gray')
title('true subspace')
% ylabel('truth')

% incoherent
i=i+1; subplot(nrows,ncols,i), cla
delsort=sort(P.del(:),'descend');
delsparse=P.del;
delsparse(delsparse<delsort(1000))=0;
imagesc(delsparse);
title('incoherent inferred subspace')


% max degree
i=i+1; subplot(nrows,ncols,i), cla
ind = get_max_edges(P.del,G);
delmax=nan(G.n);
delmax(ind)=P.del(ind);
imagesc(delmax);
title('max degree inferred subspace')


if save_stuff
    wh=[2 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'subspace_image_sim'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end
