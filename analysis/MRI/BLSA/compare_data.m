%% get BLSA32 data

% clear, clc
rootdir = '/Users/jovo/Research/data/MRI/BLSA/';
subdir  = '2010_04_14/';
files   = dir([rootdir subdir]);  files(1:5)=[];
csvfile  = importdata([rootdir subdir 'BLSA32.txt']);
s       = length(files)-1;
n       = 70;
AsX2    = nan(n,n,s);
classX2 = nan(1,s);
idsX2   = char(zeros(s,4));

for i=1:length(files)-1
    A= importdata([rootdir subdir files(i).name]);
    A(isnan(A))=0;
    AsX2(:,:,i)=A(2:n+1,2:n+1);
    idsX2(i,:)=files(i).name(16:19);
    if strcmpi(idsX2(i,3),'_')
        idsX2(i,3:4)=['0' idsX2(i,4)];
    end
    
    for k=1:length(csvfile)
        if strcmpi(idsX2(i,:),csvfile{k}(1:4))
            if strcmp(csvfile{k}(12),'m')
                classX2(i)=1; 
            else
                classX2(i)=0; 
            end
        end
    end
end

%% get BLSAX2 data

% clear, clc
rootdir = '/Users/jovo/Research/data/MRI/BLSA/will/';
subdir  = '2010_08_16/';
files   = dir([rootdir subdir]);  files(1:5)=[];
csvfile  = importdata([rootdir subdir 'BLSA32.txt']);
s       = length(files)-1;
n       = 70;
AsX2    = nan(n,n,s);
classX2 = nan(1,s);
idsX2   = char(zeros(s,4));

for i=1:length(files)-1
    A= importdata([rootdir subdir files(i).name]);
    A(isnan(A))=0;
    AsX2(:,:,i)=A(2:n+1,2:n+1);
    idsX2(i,:)=files(i).name(16:19);
    if strcmpi(idsX2(i,3),'_')
        idsX2(i,3:4)=['0' idsX2(i,4)];
    end
    
    for k=1:length(csvfile)
        if strcmpi(idsX2(i,:),csvfile{k}(1:4))
            if strcmp(csvfile{k}(12),'m')
                classX2(i)=1; 
            else
                classX2(i)=0; 
            end
        end
    end
end

%% get BLSAX2 data

% clear, clc
rootdir = '/Users/jovo/Research/data/MRI/BLSA/will/';
subdir  = '2010_09_16/';
files   = dir([rootdir subdir]);  files(1:5)=[];
csvfile  = importdata([rootdir subdir 'BLSA32.txt']);
s       = length(files)-1;
n       = 70;
AsX2    = nan(n,n,s);
classX2 = nan(1,s);
idsX2   = char(zeros(s,4));

for i=1:length(files)-1
    A= importdata([rootdir subdir files(i).name]);
    A(isnan(A))=0;
    AsX2(:,:,i)=A(2:n+1,2:n+1);
    idsX2(i,:)=files(i).name(16:19);
    if strcmpi(idsX2(i,3),'_')
        idsX2(i,3:4)=['0' idsX2(i,4)];
    end
    
    for k=1:length(csvfile)
        if strcmpi(idsX2(i,:),csvfile{k}(1:4))
            if strcmp(csvfile{k}(12),'m')
                classX2(i)=1; 
            else
                classX2(i)=0; 
            end
        end
    end
end


%%
for i=1:length(labels)
    y=labels(i).gender;
    if ~isempty(y); classX2(i)=y; end
end

%%

k=0;
num_eqs=0;
for i=1:S
    for j=1:s
        if strcmpi(ids79(i,:),idsX2(j,:))
            disp([class79(i) classX2(j)])
            if class79(i)==classX2(j), num_eqs=num_eqs+1; end
            k=k+1;
            %             figure(1)
            %             subplot(131), imagesc(As79(:,:,i)), colorbar
            %             subplot(1X2), imagesc(AsX2(:,:,j)), colorbar
            %             subplot(133), imagesc(As79(:,:,i)-AsX2(:,:,j)), colorbar
            %             keyboard
        end
    end
end

%% classify

subspace.naive_bayes=find(sum(As,3));


constants79 = get_constants(As79,class79);     % get constants to ease classification code
[Lhat79 params79 subspace79] = run_naive_bayes(As79,constants79,As79,constants79,subspace);
disp(Lhat79)

%%
constantsX2 = get_constants(AsX2,classX2);     % get constants to ease classification code
[LhatX2 paramsX2 subspaceX2] = run_naive_bayes(AsX2,constantsX2,AsX2,constantsX2,subspace);
disp(LhatX2)

