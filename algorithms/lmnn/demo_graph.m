% kidney-egg graph classification problem (one of the simplest imaginable)
clear, clc

n=10;    % number vertices
m=3;    % number of signal vertices 
p=0.4;  % probability of connection
q=0.1;  % probability of connection between signal edges in class 1
s=500; % number of training samples

E0 = p*ones(n);     % prob of connection for class 0
E1 = p*ones(n);     % prob of connection for class 1 in non signal edges
E1(1:m,1:m)=q;      % prob of connection for class 1 in signal edges

% training data
A0      = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1      = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1
siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
xTr(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
xTr(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
yTr=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

% test data
A0      = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);         % class 0 training samples
A1      = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);         % class 1
siz0=size(A0);                                          % size of array storing class 0 adjacency matrices
siz1=size(A1);                                          % size of array storing class 1 adjacency matrices
xTe(:,:,1:siz0(3))=A0;                   % create adjacency_matrices to store all sampled graphs
xTe(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;   % add class 1 samples
yTe=[zeros(1,siz0(3)) ones(1,siz1(3))];        % vector of class labels

% do embeddings
num_dim=n^2;    % number of dimensions upon embedding adjacency matrices within a vector space
xTr=double(reshape(xTr,num_dim,s));
xTe=double(reshape(xTe,num_dim,s));

%%
[L,Det]=lmnn(xTr,yTr,'quiet',1);
enerr=energyclassify(L,xTr,yTr,xTe,yTe,3);
knnerrL=knnclassify(L,xTr,yTr,xTe,yTe,3);
knnerrI=knnclassify(eye(size(L)),xTr,yTr,xTe,yTe,3);

clc;
fprintf('graph data set:\n');
fprintf('3-NN Euclidean training error: %2.2f\n',knnerrI(1)*100);
fprintf('3-NN Euclidean testing error: %2.2f\n',knnerrI(2)*100);
fprintf('3-NN Malhalanobis training error: %2.2f\n',knnerrL(1)*100);
fprintf('3-NN Malhalanobis testing error: %2.2f\n',knnerrL(2)*100);
fprintf('Energy classification error: %2.2f\n',enerr*100);
fprintf('Training time: %2.2fs\n\n\n',Det.time);

