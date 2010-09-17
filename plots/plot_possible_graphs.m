figure(1), clf
n=1:50;
y=0*n;
for i=n
    y(i)=2^choose(i,2);
end
h=semilogy(n,y,'k','linewidth',2);
fs=12;
xlabel('number of vertices','fontsize',fs)
ylabel('number of distinct simple graphs','fontsize',fs)

%%
wh=[3 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=['num_of_graphs'];
print('-dpdf',figname)