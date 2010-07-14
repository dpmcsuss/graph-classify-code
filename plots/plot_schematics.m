%% number of simple graphs as a function of number of nodes
figure(1), clf
n=1:150;
y=0*n;
for i=n
    y(i)=choose(i,2);
end
h=plot(n,y,'k','linewidth',2);
fs=12;
xlabel('number of vertices','fontsize',fs)
ylabel('number of distinct simple graphs','fontsize',fs)

wh=[2.5 4];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname='~/Research/necog/figs/misc/num_of_graphs';
print('-dpdf',figname)


%% plot xor
figure(2), clf, clc
s=50;
x=0.35*randn(8,s);
x=x+repmat([0 0 1 1 0 1 1 0]',1,s);
z=rand(2,1)+[0.5; 0.5];

xTr = [x(1:2,:) x(3:4,:) x(5:6,:) x(7:8,:)];
yTr=[zeros(1,s*2) ones(1,s*2)];
[L,Det]=lmnn(xTr,yTr,'quiet',1);
xnew=L*xTr;

subplot(121), cla, hold all
scatter(xTr(1,1:s*2),xTr(2,1:s*2),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
scatter(xTr(1,s*2+1:end),xTr(2,s*2+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);

scatter(z(1),z(2),'filled','markerfacecolor','k')
circle(z,0.18,500,'k-');
axis([-.5 1.5 -.5 1.5])

subplot(122), cla, hold all
scatter(xnew(1,1:s*2),xnew(2,1:s*2),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
scatter(xnew(1,s*2+1:end),xnew(2,s*2+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);

znew=L*z;
scatter(znew(1),znew(2),'filled','markerfacecolor','k')
circle(znew,0.18,500,'k-')
axis('tight')
% axis([-.5 1.5 -.5 1.5])

