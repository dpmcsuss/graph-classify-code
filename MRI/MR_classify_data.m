clear, clc, clf, %close all

datadir = '~/Research/data/MRI/BLSA/will/';
figdir  = '~/Research/figs/MRI/BLSA/';
fname   = 'BLSA42';
load([datadir fname])

save_stuff = 0;

%% prepare data
for i=1:length(labels)
    A=As(:,:,i);
    Y = 0.2;                    % binarize 
    %Y = quantile(A(:),0.95);
    A(A<=Y)=0;
    As(:,:,i)=A;
    ys(i)=labels(i).gender;     % get Y from labels
end

G = get_constants(As,ys);

if save_stuff, save([datadir fname '_params']); end

%% in sample

Z.case  = 'in';
Z.nb    = 1;
Z.max   = 1; Z.Nmax = G.n/10;
Z.inc   = 1; Z.Ninc = G.n;

[Lhatin ind Pin yhatin] = graph_classify_ind_edge(As,G,Z);

%% out-of-sample

nmin=min(G.n0,G.n1);
nmax=max(G.n0,G.n1);

clear Z
Z.nb    = 1;
Z.inc   = 1; Z.Ninc = G.n;
Z.max   = 1; Z.Nmax = G.n/10;
Z.case  = 'lso';

Ntrain=[2 4 6 8 10 12 14];
Nresamples=100;

Lhats.nb=zeros(length(Ntrain),Nresamples);
Lhats.inc=zeros(length(Ntrain),Nresamples);
Lhats.max=zeros(length(Ntrain),Nresamples);

for i=1:length(Ntrain)

    for j=1:Nresamples
        
        ind0  = randperm(nmin);
        ind1  = randperm(nmin);
        
        y0trn = G.y0(ind0(1:Ntrain(i)));
        y0tst = G.y0(ind0(Ntrain(i)+1:end)); 
        
        y1trn = G.y1(ind1(1:Ntrain(i)));
        y1tst = G.y1(ind1(Ntrain(i)+1:end));
        
        Atrn = As(:,:,[y0trn y1trn]);
        Atst = As(:,:,[y0tst y1tst]);
        
        ytrn = [zeros(1,Ntrain(i)) ones(1,Ntrain(i))];
        ytst = [zeros(1,length(y1tst)) ones(1,length(y1tst))];
        
        Gtrn = get_constants(Atrn,ytrn);
        Gtst = get_constants(Atst,ytst);
                
        Lhat = graph_classify_ie(Atrn,Gtrn,Z,Atst,Gtst);

        Lhats.nb(i,j)=Lhat.nb;
        Lhats.inc(i,j)=Lhat.inc;
        Lhats.max(i,j)=Lhat.max;

    end
end

if save_stuff, save([datadir fname]); end

%% plot params

figure(1), clf
nrows=1;
ncols=3;

P = get_params(As,G);

emax=max(max([P.E0(:) P.E1(:) abs(P.E0(:)-P.E1(:))]));
subplot(nrows,ncols,1), cla
image(60*P.E0/emax)
colormap('gray')
title('average female')
xlabel('vertices')
ylabel('vertices')

subplot(nrows,ncols,2), cla
image(60*P.E1/emax)
title('average male')

subplot(nrows,ncols,3), cla
image(60*abs(P.E0-P.E1)/emax)
colormap('gray')
title('difference')

if save_stuff
    wh=[6 1.8];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir fname '_params'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% plot misclassification
figure(888), clf, hold all

ls1='-'; ls2='--'; ls3=':'; ls4='.';
color1=0*[1 1 1]; color2=0.75*[1 1 1]; color3=0.7*[1 1 1]; color4=0.25*[1 1 1];
lw=2;
fs=12;
xticks=Ntrain;
clear Ncorrect ind_IE_mat ind_max_mat

h(1)=errorbar(Ntrain,mean(Lhats.nb,2),var(Lhats.nb,[],2),'color',color1,'linestyle',ls1,'linewidth',lw);
h(2)=errorbar(Ntrain+0.2,mean(Lhats.inc,2),var(Lhats.inc,[],2),'color',color2,'linestyle',ls2,'linewidth',lw);
h(3)=errorbar(Ntrain+0.4,mean(Lhats.max,2),var(Lhats.max,[],2),'color',color3,'linestyle',ls3,'linewidth',lw);
% h(4)=plot([0 length(labels)],[0.5 0.5],'--','color',0.5*[1 1 1]);
axis([0 max(Ntrain)+1 0 0.5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples per class','fontsize',fs)
legend(gca,'unconstrained','edge constrained','vertex constrained','Location','SouthWest')
set(gca,'XTick',Ntrain,'XTickLabel',Ntrain)


if save_stuff
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir fname '_rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% plot example recovery image

load('~/Research/necog/data/BLSA/will/BLSA32')
figure(11), clf
nrows=3;
ncols=1;
i=0;

G = get_constants(As,ys);
P = get_params(As,G);

% true signal
i=i+1; subplot(nrows,ncols,i), cla
imagesc(P.del)
colormap('gray')
title('unconstrained')
% ylabel('truth')

% incoherent
i=i+1; subplot(nrows,ncols,i), cla
delsort=sort(P.del(:),'descend');
delsparse=P.del;
delsparse(delsparse<delsort(G.n))=0;
imagesc(delsparse);
title('edge constrained subspace')

% max degree
i=i+1; subplot(nrows,ncols,i), cla
G.Nmax=G.n/10;
ind = get_max_edges(P.del);
delmax=zeros(G.n);
delmax(ind)=P.del(ind);
imagesc(delmax);
title('vertex constrained subspace')


if save_stuff
    wh=[2 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir fname '_subspaces'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end
