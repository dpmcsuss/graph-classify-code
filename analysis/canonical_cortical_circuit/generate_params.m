clear, clc

[num,txt,raw]=xlsread('/Users/jovo/Research/data/canonical_cortical_circuit/matrix_of_connectivity.csv');

AdjMat=zeros(8);

n=length(txt);

L1=3-1;
L2e=[4:5 17]-1;
L2i=[9:16 18]-1;
L4e=[19:22]-1;
L4i=[23:33]-1;
L5e=[34:42]-1;
L5i=[43:48]-1;
L6=[49:55]-1;

for i=2:n
    if i==L1,           I=1;    % layer 1 (inhibitory)
    elseif any(i==L2e), I=2;    % layer 2/3 excitatory
    elseif any(i==L2i), I=3;    % layer 2/3 inhibitory
    elseif any(i==L4e), I=4;    % layer 4 excitatory
    elseif any(i==L4i), I=5;    % layer 4 inhibitory
    elseif any(i==L5e), I=6;    % layer 5 excitatory
    elseif any(i==L5i), I=7;    % layer 5 inhibitory
    elseif any(i==L6),  I=8;    % layer 6 (excitatory)
    end
        
    for j=2:n
        if j==L1,           J=1;    % layer 1 (inhibitory)
        elseif any(j==L2e), J=2;    % layer 2/3 excitatory
        elseif any(j==L2i), J=3;    % layer 2/3 inhibitory
        elseif any(j==L4e), J=4;    % layer 4 excitatory
        elseif any(j==L4i), J=5;    % layer 4 inhibitory
        elseif any(j==L5e), J=6;    % layer 5 excitatory
        elseif any(j==L5i), J=7;    % layer 5 inhibitory
        elseif any(j==L6),  J=8;    % layer 6 (excitatory)
        end
    
        AdjMat(I,J)=AdjMat(I,J)+raw{i,j};
        
    end
end


% block structure
B=abs(AdjMat);
B=B/(17*mean(B(:)));

w(1) = 1.5;         % layer 1 (inh)
w(2) = 26;          % layer 2/3 exc
w(3) = 7.3;         % layer 2/3 inh
w(4) = 9.2+9.2+9.2; % layer 4 e
w(5) = 6.9;         % layer 4 i
w(6) = 4.8+1.3;     % layer 5 e
w(7) = 0.6+0.8;     % layer 5 i
w(8) = 13.6+4.5;    % layer 6 e
% w(9) = 4;           % layer 6 i

w=w/sum(w);

w*B*w'


save('/Users/jovo/Research/data/canonical_cortical_circuit/params.mat','B','w')
