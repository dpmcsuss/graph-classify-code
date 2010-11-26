function [d_ind deg_ind] = plot_subspaces(del,subspace,n,alg)

fsx=8;
fs=10;
lw=2;

figure(2), 
subplot(2,3,1)
imagesc(del);
colormap('gray')
title('absolute difference','fontsize',fs)
% xlabel('vertices','fontsize',fs)
% ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'fontsize',fsx)

subplot(2,3,2)
[d_val d_ind] = sort(del(:),'descend');
d_ind=d_ind(1:choose(n,2));
plot(d_val(1:choose(n,2)),'k','linewidth',lw)
axis('tight')
xlabel('sorted edges','fontsize',fsx)
ylabel('$\hat{\delta}_{ij}$','fontsize',fsx,'interp','latex')
title('edge differences','fontsize',fs)
set(gca,'XTick',[0 1000 2000],'fontsize',fsx)

subplot(2,3,3)
[deg_val deg_ind] = sort(sum(del),'descend');
plot(deg_val,'k','linewidth',lw)
title('degree differences','fontsize',fs)
xlabel('sorted vertices','fontsize',fsx)
ylabel('$\hat{\Delta}_{ij}$','fontsize',fsx,'interp','latex')
axis('tight')
set(gca,'fontsize',fsx)


subplot(2,3,4)
sub.all=zeros(n);
sub.all(subspace.naive_bayes)=1;
imagesc(sub.all)
colormap('gray')
title('naive bayes','fontsize',fs)
xlabel('vertices','fontsize',fsx)
ylabel('vertices','fontsize',fsx)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'fontsize',fsx)

subplot(2,3,5)
sub.inc=zeros(n);
sub.inc(d_ind(1:n))=1;
imagesc(sub.inc)
colormap('gray')
title('incoherent','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
set(gca,'fontsize',fsx)

subplot(2,3,6)
sub.star=zeros(n);
deg_ind = get_max_edges(del,8);
sub.star(deg_ind)=1;
sub.star=triu(sub.star,1);
imagesc(sub.star)
colormap('gray')
title('coherent','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])


if alg.save
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_subspaces'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
