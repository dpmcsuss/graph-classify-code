clear, clc, clf
rand('seed',12346);

load('~/Research/necog/data/c_elegans/c_elegans_chemical_connectome')

n       = length(A);        % |V|=n
A(isnan(A))=0;          % remove NaN's
A       = A+A';         % fill in lower triangle (this is an undirected graph)

s       = 200;          % # of samples
eps     = 0.05;         % noise parameter "eps"

% make egg
thesem  = [69,80,82,94,110,127,129,133,134,138];% which vertices are the signal containing vertices
m       = length(thesem);            % size of egg
sig     = 2;            % egg parameter "sig"
egg     = rand(m)*sig-sig/2;

% make E0
E0      = A + eps*rand(n);
E0(E0<=0)   = eps;

% make E1 = E0 + egg
E1          = E0;
E1(thesem,thesem) = E1(thesem,thesem) + egg;
E1(E1<=0)   = eps;

Bs  = zeros(n,n,s);
Bs(:,:,1:s/2)   = poissrnd(repmat(E0,[1 1 s/2]));
Bs(:,:,s/2+1:s) = poissrnd(repmat(E1,[1 1 s/2]));


ind_mat=zeros(n);
ind_mat(thesem,thesem)=1;
trueind=find(ind_mat);

datadir = '~/Research/necog/data/c_elegans/';
figdir  = '~/Research/necog/figs/c_elegans/';
fname   = 'egg_sim';
save_stuff = 1;

if save_stuff, save([datadir fname]); end

%% plot sim

figure(1), clf
nrows=1;
ncols=3;
xticks=[70 140 210 279];

emax=max(max([E0(:) E1(:) abs(E0(:)-E1(:))]));
subplot(nrows,ncols,1), cla
image(60*E0/emax)
colormap('gray')
title('E_0')
xlabel('pre-synaptic neuron')
ylabel('post-synaptic neuron')
set(gca,'XTick',xticks,'YTick',xticks)

subplot(nrows,ncols,2), cla
image(60*E1/emax)
title('E_1')
set(gca,'XTick',xticks,'YTick',xticks)

subplot(nrows,ncols,3), cla
image(60*abs(E0-E1)/emax)
colormap('gray')
title('\Delta')
set(gca,'XTick',xticks,'YTick',xticks)

if save_stuff
    wh=[6 1.8];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'P0_P1_Del'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% in sample

ys = [zeros(1,s/2) ones(1,s/2)];
G = get_constants(Bs,ys);

G.nb    = [];
G.inc   = []; G.Ninc=m^2;
G.max   = []; G.Nmax=m;

[Lhatin ind P yhatin] = graph_classify(Bs,G);

%% out-of-sample

k=0;
is=25:25:s/2; %10:10:50; %[10 25 50 75 100];
for i=is

    BBs=zeros(n,n,2*i);
    BBs(:,:,1:i)=Bs(:,:,1:i);
    BBs(:,:,i+1:2*i)=Bs(:,:,s/2+1:s/2+i);

    ys = [zeros(1,i) ones(1,i)];
    G = get_constants(BBs,ys);

    G.nb    = [];
    G.inc   = []; G.Ninc=m^2;
    G.max   = []; G.Nmax=m;
    G.loo   = [];
    G.tru   = trueind;

    [Lhat ind] = graph_classify(BBs,G);

    k=k+1;

    Lhats.nb(k)=Lhat.nb;
    Lhats.inc(k)=Lhat.inc;
    Lhats.max(k)=Lhat.max;
    Lhats.tru(k)=Lhat.tru;

    inds{k}  = ind;

end

if save_stuff, save([datadir fname]); end

%% plot example recovery image

figure(11), clf
nrows=3;
ncols=2;
i=0;

% true signal
i=i+1; subplot(nrows,ncols,i), cla
imagesc(ind_mat)
colormap('gray')
title('true subspace')
% ylabel('truth')

i=i+1; subplot(nrows,ncols,i), cla
imagesc(ind_mat)
axis([thesem(1)-.5 thesem(end)+.5 thesem(1)-.5 thesem(end)+.5])
title('zoomed subspaces')

% naive bayes
i=i+1; subplot(nrows,ncols,i), cla
ind_IE_mat = zeros(G.n);
ind_IE_mat(inds{4}(1).inc) = 1;
imagesc(ind_IE_mat);
title('incoherent inferred subspace')

i=i+1; subplot(nrows,ncols,i), cla
imagesc(ind_IE_mat);
axis([thesem(1)-.5 thesem(end)+.5 thesem(1)-.5 thesem(end)+.5])

% max degree
i=i+1; subplot(nrows,ncols,i), cla
ind_max_mat = zeros(G.n);
ind_max_mat(inds{4}(1).max) = 1;
imagesc(ind_max_mat);
title('max degree inferred subspace')


i=i+1; subplot(nrows,ncols,i), cla
imagesc(ind_max_mat);
axis([thesem(1)-.5 thesem(end)+.5 thesem(1)-.5 thesem(end)+.5])

if save_stuff
    wh=[4 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'subspace_image_sim'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

%% recovered signals

figure(4), clf
ms = 12;
gray = .75*[1 1 1];

sigdims = trueind;
i       = max(is);
BBs     = zeros(n,n,2*i);
BBs(:,:,1:i)    = Bs(:,:,1:i);
BBs(:,:,i+1:2*i)= Bs(:,:,s/2+1:s/2+i);
ys      = [zeros(1,i) ones(1,i)];
G       = get_constants(BBs,ys);
P       = get_params(BBs,G,i);
delhat  = abs(P.E0-P.E1);

subplot(121)
nedges=100;
delsort = sort(delhat(:),'descend');
plot(delsort(1:nedges),'.k','markersize',ms)
hold all
delsort2=nan(size(delhat));
for i=1:n^2, if any(delsort(i)==delhat(sigdims)), delsort2(i)=delsort(i); end; end
plot(delsort2(1:nedges),'.','color',gray,'markersize',ms)
hold off
title(['incoherent edge recovery'])
ylab=ylabel('$\hat{\delta}_{ij}$','interp','latex');
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('sorted edge index')
axis('tight')
Nedges_correct = sum(sum(ind_IE_mat(thesem,thesem)));

text(10,max(delhat(:))*0.95,['% correct edges = ' num2str(round(100*Nedges_correct/m^2)/100)])


% max degree
G.Nmax = m;
ind_max = get_max_edges(delhat,G);
deg = sum(delhat,1) + sum(delhat,2)';
degtrue=nan(n,1);
degtrue(thesem)=deg(thesem);
[degsort IX] = sort(deg,'descend');
deg_correct = find(deg(thesem)>=degsort(m));
Ndeg_correct = length(deg_correct);

subplot(122)
hold all;
[foo bar]=sort(deg,'descend');
plot(deg(bar),'.k','markersize',ms)
plot(degtrue(bar),'.','color',gray,'markersize',ms)
axis('tight')
ylab=ylabel('estimated degree');
xlabel('sorted vertex index')
title('max degree vertex recovery')
box on

text(50,max(deg)*.95,['% correct vertices = ' num2str(round(100*Ndeg_correct/m)/100)])

wh=[6 2];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'subspace_plots_sim'];
print('-dpdf',figname)
saveas(gcf,figname)


%% plot misclassification and edge detection rates
figure(888), clf, hold all

ls1='-'; ls2='--'; ls3=':'; ls4='.';
color1=0*[1 1 1]; color2=0.75*[1 1 1]; color3=0.7*[1 1 1]; color4=0.25*[1 1 1];
lw=2;
fs=12;
xticks=is;
clear Ncorrect ind_IE_mat ind_max_mat

subplot(121), hold all
h(1)=plot(is,Lhats.nb,'color',color1,'linestyle',ls1,'linewidth',lw);
h(2)=plot(is,Lhats.inc,'color',color2,'linestyle',ls2,'linewidth',lw);
h(3)=plot(is,Lhats.max,'color',color3,'linestyle',ls3,'linewidth',lw);
h(4)=plot(is,Lhats.tru,'color',color4,'linestyle',ls4,'linewidth',lw);
axis([0 max(is) 0 max(Lhats.nb)])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)
legend(h,'nb','inc','max','true','location','best')
set(gca,'XTick',is,'XTickLabel',2*is)


subplot(122), hold all
for k=1:length(is)
    jmax=length(inds{k});
    for j=1:jmax
        Zmat = zeros(G.n);        
        Zmat(inds{k}(j).inc) = 1;
        Ncorrect.inc(j) = sum(sum(Zmat(trueind)));

        Zmat = zeros(G.n);
        Zmat(inds{k}(j).max) = 1;
        Ncorrect.max(j) = sum(sum(Zmat(trueind)));
    end
    mean_correct.inc(k) = mean(Ncorrect.inc)/m^2;
    std_correct.inc(k) = std(Ncorrect.inc)/m^2;

    mean_correct.max(k) = mean(Ncorrect.max)/m^2;
    std_correct.max(k) = std(Ncorrect.max)/m^2;
end

errorbar(mean_correct.inc,std_correct.inc,'color',color2,'linestyle',ls2,'linewidth',lw)
errorbar(mean_correct.max,std_correct.max,'color',color3,'linestyle',ls3,'linewidth',lw)

axis([0 k 0 1])
set(gca,'XTick',1:k,'XTickLabel',2*is)
ylabel('edge detection rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)

if save_stuff
    wh=[5 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end

