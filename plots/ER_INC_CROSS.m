n=10;
m=8;
p=0.2;
q=0.5;



M_inc=randperm(n^2);
M_inc=M_inc(1:m);

M_cross=randperm(n);
M_cross=M_cross(1:m);


ER=p*ones(n);

INC=p*ones(n);
INC(M_inc)=q;

CROSS=p*ones(n);
CROSS(5,M_cross(1:m/2))=q;
CROSS(M_cross(m/2+1:m),5)=q;

figure(1), clf
if nargin==3, nrows=2; else nrows=1; end
ncols=3;
fs=16;


emax=q;
subplot(131), cla
image(60*ER/emax)
colormap('gray')
title('Erdos-Renyi','fontsize',fs)
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

subplot(132), cla
image(60*INC/emax)
title('Incoherent','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])

subplot(133), cla
image(60*CROSS/emax)
colormap('gray')
title('Star_1','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])


%%
wh=[7.5 3.5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname='../../../figs/misc/ER_INC_STAR1';
print('-dpdf',figname)
saveas(gcf,figname)


