
% load('data/bal2.mat');

s=1000;
x=0.2*randn(8,s);
x=x+repmat([0 0 1 1 0 1 1 0]',1,s);

xTr = [x(1:2,1:s/2) x(3:4,1:s/2) x(5:6,1:s/2) x(7:8,1:s/2)];
yTr=[zeros(1,s) ones(1,s)];

xTe = [x(1:2,1+s/2:end) x(3:4,1+s/2:end) x(5:6,1+s/2:end) x(7:8,1+s/2:end)];
yTe=[zeros(1,s) ones(1,s)];

%%
[L,Det]=lmnn(xTr,yTr,'quiet',1);
enerr=energyclassify(L,xTr,yTr,xTe,yTe,3);
knnerrL=knnclassify(L,xTr,yTr,xTe,yTe,3);
knnerrI=knnclassify(eye(size(L)),xTr,yTr,xTe,yTe,3);

clc;
fprintf('Bal data set:\n');
fprintf('3-NN Euclidean training error: %2.2f\n',knnerrI(1)*100);
fprintf('3-NN Euclidean testing error: %2.2f\n',knnerrI(2)*100);
fprintf('3-NN Malhalanobis training error: %2.2f\n',knnerrL(1)*100);
fprintf('3-NN Malhalanobis testing error: %2.2f\n',knnerrL(2)*100);
fprintf('Energy classification error: %2.2f\n',enerr*100);
fprintf('Training time: %2.2fs\n (As a reality check: My desktop needs 20.53s)\n\n',Det.time);

xnew=L*xTr;

%% 
subplot(121), cla, hold all
scatter(xTr(1,1:s),xTr(2,1:s),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
scatter(xTr(1,s+1:end),xTr(2,s+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);
axis('tight')

subplot(122), cla, hold all
scatter(xnew(1,1:s),xnew(2,1:s),'marker','+','markerfacecolor','r','markeredgecolor','r','linewidth',2);
scatter(xnew(1,s+1:end),xnew(2,s+1:end),'marker','.','markerfacecolor','b','markeredgecolor','b','sizedata',300);
axis('tight')
