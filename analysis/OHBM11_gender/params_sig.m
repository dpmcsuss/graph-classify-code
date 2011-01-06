
clc
figure(1), clf,
fs=10;

load('/Users/jovo/Research/data/MRI/BLSA/results/edge_significance.mat')


etc.save=0;
etc.figdir = '/Users/jovo/Research/figs/MRI/BLSA/';
etc.fname = 'BLSA50';


params = get_params_mle(As,constants);

E0 = tril(params.E0+params.E0');
E1 = tril(params.E1+params.E1');

emax=max(max([params.E0(:) params.E1(:) abs(params.E0(:)-params.E1(:))]));

subplot(221)
image(64*(1-E0)/emax)
colormap('gray')
title('female','fontsize',fs)
xlabel('gyral indices','fontsize',fs)
ylabel('gyral indices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'XTick',[n/2 n],'YTick',[n/2 n])
box off

subplot(222)
image(64*(1-E1)/emax)
title('male','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'XTick',[n/2 n],'YTick',[n/2 n])
box off

subplot(224)
sym=tril((1-d_pval)+(1-d_pval'))+triu(ones(n));
imagesc(sym)
title('significance of difference','fontsize',fs)
set(gca,'fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'XTick',[n/2 n],'YTick',[n/2 n])
box off
axis square

A = imread('/Users/jovo/Research/figs/MRI/misc/Adjacency_brain_figure.jpg');
subplot(223)
imagesc(A);
set(gca,'XTick',[],'YTick',[])


etc.save=1;
if etc.save
    wh=[4 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[etc.figdir etc.fname '_params'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
