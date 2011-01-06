clear, clc

rootdir    = '/Users/jovo/Research/data/MRI/BLSA/';
datadir{1} = '2010_04_14/';
datadir{2} = '2010_05_27/';
datadir{3} = '2010_08_16/';
datadir{4} = '2010_09_16/';
datadir{5} = '2010_11_07/';
datadir{6} = '2010_11_07_dups/';


labels = importdata([rootdir 'blsa_subs_exactAges.csv'],',');


for i=1:length(datadir)
    files{i} = dir([rootdir datadir{i}]);
    files{i}(1:2)=[];
    if strcmpi(files{i}(1).name,'.DS_Store'), files{i}(1)=[]; end
end

%% 4/14 data

j=1;
k=0;
for i=1:length(files{j})
    k=k+1;
    exp(k).fname = [datadir{j} files{j}(i).name];
    A = importdata([rootdir exp(k).fname]);
    exp(k).Adj = A(2:end,2:end);
    exp(k).id  = [files{j}(i).name(16:17) files{j}(i).name(13:14)];

    for l=2:length(labels)
        if strcmpi(labels{l}(1:4),exp(k).id)
            if strcmpi(labels{l}(end-3),'M') || strcmpi(labels{l}(end-2),'M')
                exp(k).sex = 1;
            elseif strcmpi(labels{l}(end-3),'F') || strcmpi(labels{l}(end-2),'F')
                exp(k).sex = 0;
            else
                keyboard
            end
        end
    end
end

%% 5/27, 8/16, and 9/16 data

for j=2:4
    for i=1:length(files{j})
        k=k+1;
        exp(k).fname = [datadir{j} files{j}(i).name];
        A = importdata([rootdir exp(k).fname]);
        exp(k).Adj = A(2:end,2:end);
        exp(k).id  = files{j}(i).name(16:19);

        for l=2:length(labels)
            if strcmpi(labels{l}(1:4),exp(k).id)
                if strcmpi(labels{l}(end-3),'M') || strcmpi(labels{l}(end-2),'M')
                    exp(k).sex = 1;
                elseif strcmpi(labels{l}(end-3),'F') || strcmpi(labels{l}(end-2),'F')
                    exp(k).sex = 0;
                else
                    keyboard
                end
            end
        end
    end
end

%% 11/7 good data

j=5;
for i=1:length(files{j})
    k=k+1;
    exp(k).fname = [datadir{j} files{j}(i).name];
    A = importdata([rootdir exp(k).fname]);
    exp(k).Adj = A(2:end,2:end);
    exp(k).id  = files{j}(i).name(24:27);

    for l=2:length(labels)
        if strcmpi(labels{l}(1:4),exp(k).id)
            if strcmpi(labels{l}(end-3),'M') || strcmpi(labels{l}(end-2),'M')
                exp(k).sex = 1;
            elseif strcmpi(labels{l}(end-3),'F') || strcmpi(labels{l}(end-2),'F')
                exp(k).sex = 0;
            else
                keyboard
            end
        end
    end
end


%% 11/7 duplicate data

j=6;
for i=1:length(files{j})
    k=k+1;
    exp(k).fname = [datadir{j} files{j}(i).name];
    A = importdata([rootdir exp(k).fname]);
    exp(k).Adj = A(2:end,2:end);
    exp(k).id  = files{j}(i).name(29:32);

    for l=2:length(labels)
        if strcmpi(labels{l}(1:4),exp(k).id)
            if strcmpi(labels{l}(end-3),'M') || strcmpi(labels{l}(end-2),'M')
                exp(k).sex = 1;
            elseif strcmpi(labels{l}(end-3),'F') || strcmpi(labels{l}(end-2),'F')
                exp(k).sex = 0;
            else
                keyboard
            end
        end
    end
end


%% compare duplicates

clear A comp keepers
cont=0;

a=ones(70);
a=triu(a,+1);
nb=find(a);
keepers(1)=exp(1);
ll=0;
l=0;
numwithzero=0;
for k=1:length(exp)
    clc
    disp(['k=' num2str(k)])
    for kkk=1:k-1
        if strcmpi(exp(k).id,exp(kkk).id)
            %             exp(k)
            %             exp(kkk)
            %             disp('already did comp')
            %             keyboard
            cont=1;
            break
        end
    end
    if cont, cont=0; continue, end
    l=l+1; disp(['l=' num2str(l)])
    j=1;
    comp{k}(j)=exp(k);
    exp(k)
    for kk=k+1:length(exp)
        if strcmpi(exp(k).id,exp(kk).id)
            j=j+1;
            comp{k}(j)=exp(kk);
            exp(kk)
        end
    end
    AA=comp{k}(1).Adj;
    AA(isnan(AA))=0;
    for jj=1:j
        BB=comp{k}(jj).Adj;
        BB(isnan(BB))=0;
        subplot(1,j,jj), imagesc(BB);
        err = num2str(norm(AA-BB));
        nzeros = num2str(length(find(comp{k}(jj).Adj(nb)==0)));
        if (jj==1 && str2num(nzeros)>0), numwithzero=l; end
        nl2 = num2str(length(find(comp{k}(jj).Adj(nb)<0.2)));
        title(['err=' err ', #0=' nzeros ', <.2=' nl2])
    end
    keyboard
end

%% keep only unique keepers

clear uniq_keeps
uniq_keeps(1)=keepers(1);
kkk=0; cont=0;
for k=2:length(keepers)
    kkk=kkk+1;
    uniq_keeps(kkk)=keepers(k)
    for kk=2:k-1
        if strcmpi(keepers(k).id(1:2),keepers(kk).id(1:2)),
            kkk=kkk-1;
            cont=1;
            break
        end
    end
    if cont, cont=0; continue, end    
end


save([rootdir 'results/validated_data'])
save([rootdir 'results/BLSA50_data'],'uniq_keeps')
