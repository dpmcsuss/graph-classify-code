%% number of simple graphs as a function of number of nodes
figure(1), clf
n=1:50;
y=0*n;
for i=n
    y(i)=choose(i,2);
end
h=semilogy(n,2.^y,'k','linewidth',2);
fs=12;
xlabel('number of vertices','fontsize',fs)
ylabel('number of distinct simple graphs','fontsize',fs)

wh=[4 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname='~/Research/necog/figs/misc/num_of_graphs';
print('-dpdf',figname)


%% plot xor
figure(2), clf, clc
s=50;
x=0.2*randn(8,s);
x=x+repmat([-.5 -.5 .5 .5 -.5 .5 .5 -.5]',1,s);

xTr = [x(1:2,:) x(3:4,:) x(5:6,:) x(7:8,:)];
yTr=[zeros(1,s*2) ones(1,s*2)];
xnew=L*xTr;

% subplot(121), 
%%
cla, 
hold all
scatter(xTr(1,1:s*2),xTr(2,1:s*2),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
scatter(xTr(1,s*2+1:end),xTr(2,s*2+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);

z=[0.2; 0];
scatter(z(1),z(2),'marker','*','markerfacecolor','k','markeredgecolor','k','sizedata',72)
circle(z,0.21,500,'k-');
axis([-1 1 -1 1])
xlabel('dimension 1')
ylabel('dimension 2')

wh=[3 3]*1.0;   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
alg.figdir  = '~/Research/necog/figs/misc/';
alg.fname   = 'knn_schem';

figname=[alg.figdir alg.fname];
print('-dpdf',figname)
saveas(gcf,figname)


% [L,Det]=lmnn(xTr,yTr,'quiet',1);
% subplot(122), cla, hold all
% scatter(xnew(1,1:s*2),xnew(2,1:s*2),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
% scatter(xnew(1,s*2+1:end),xnew(2,s*2+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);
% 
% znew=L*z;
% scatter(znew(1),znew(2),'filled','markerfacecolor','k')
% circle(znew,0.18,500,'k-')
% axis('tight')
% axis([-.5 1.5 -.5 1.5])

%% plot independe edge schematic
clear; clc

n=9;

E0=rand(n);
E0=E0+E0';
diags=1:n+1:n^2;
E0(diags)=0;

m=2;
subgraph=randperm(n^2);
subgraph=subgraph(1:m);
[I J]=ind2sub([n n],subgraph);

E1=E0;
k=0;
rnd=rand(100,1);
for i=I
    for j=J
        k=k+1;
        E1(i,j)=rnd(k);
        E1(j,i)=rnd(k);
    end
end

figure(1), clf
subplot(131), imagesc(E0), title('class 0'), xlabel('vertex number'), ylabel('vertex number'), 
set(gca,'XTick',3:3:9,'YTick',1:3:9,'YTickLabel',[9 6 3])
subplot(132), imagesc(E1), title('class 1'), set(gca,'XTicklabel',[],'YTicklabel',[])
subplot(133), imagesc(abs(E0-E1)), title('difference'), set(gca,'XTicklabel',[],'YTicklabel',[])
colormap('gray')

wh=[3 1.5]*1.0;   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
alg.figdir  = '~/Research/necog/figs/misc/';
alg.fname   = 'ie_schem';
figname=[alg.figdir alg.fname];
print('-dpdf',figname)
saveas(gcf,figname)

%% plot independe edge schematic
clear; clc

n=9;

diffuse=rand(n);
diffuse=diffuse+diffuse';
diags=1:n+1:n^2;
diffuse(diags)=0;


m=3;
subgraph=randperm(n^2);
subgraph=subgraph(1:m);
[I J]=ind2sub([n n],subgraph);

E1=zeros(9);
k=0;
rnd=rand(100,1);
for i=I
    for j=J
        k=k+1;
        E1(i,j)=diffuse(i,j);
        E1(j,i)=diffuse(j,i);
    end
end

figure(1), clf
subplot(131), 
imagesc(diffuse), 
title('dense'), 
xlabel('vertex number'), ylabel('vertex number'), 
set(gca,'XTick',3:3:9,'YTick',1:3:9,'YTickLabel',[9 6 3])

subplot(132), 
imagesc(E1), title('sparse'), set(gca,'XTicklabel',[],'YTicklabel',[])


subplot(133), 
coherent=zeros(n);
coh_graph=[1:3 10:13 19:22];
coherent(coh_graph)=diffuse(coh_graph);
imagesc(coherent)

title('structured'), set(gca,'XTicklabel',[],'YTicklabel',[])
colormap('gray')


wh=[3 1.5]*1.0;   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
alg.figdir  = '~/Research/necog/figs/misc/';
alg.fname   = 'subgraph_types';
figname=[alg.figdir alg.fname];
print('-dpdf',figname)
saveas(gcf,figname)


