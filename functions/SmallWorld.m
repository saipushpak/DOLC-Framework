%%This function is to generate undirected small-world network and plot it.
% 2013-11-05 by Xi Zhang
%Inputs:
%N: the number of nodes in network
%p: rewiring probability p, which determines the probability of switching one end of
%   each link to connect it with a randomly selected node.
%m: initialized m-regular graph
%Outputs:
%A_adj:the adjacency matrix
%A_weight: the adjacency matrix with weights for each link

function A_adj = SmallWorld(N,p,m)

% ===================== get the adjancency matrix=====================
% set up a N*N sparse matrix
matrix=sparse([],[],[],N,N,0);
for i=m+1:N-m
    for j=i-m:i+m
        matrix(i,j)=1;
    end
end
for i=1:m
    for j=1:i+m
        matrix(i,j)=1;
    end
end
for i=N-m+1:N
    for j=i-m:N
%         disp(i)
%         disp(j)
        matrix(i,j)=1;
    end
end
for i=1:m
    for j=N-m+i:N
        matrix(i,j)=1;
        matrix(j,i)=1;
    end
end
% ==========================================================

%reconnect the edge anticlockwisely from 1 to N-m-1
for i=1:N-m-1
    for j=i+1:i+m
        % choose a value randomly
        r=rand(1);
        if r<=p
            % find the nonzero
            unconect=find(matrix(i,:)==0);
            % find the number of nonzero
            M=length(unconect);
            r1=ceil(M*rand(1));
            % connect the pair
            matrix(i,unconect(r1))=1;
            matrix(unconect(r1),i)=1;
            matrix(i,j)=0;
            matrix(j,i)=0;
        end
    end
end
%reconnect the edge anticlockwisely from  N-m to N-1
for i=N-m+1:N-1
    for j=[i+1:N 1:i- N+m]
        r=rand(1);
        if r<=p
            unconect=find(matrix(i,:)==0);
            r1=ceil(length(unconect)*rand(1));
            matrix(i,unconect(r1))=1;
            matrix(unconect(r1),i)=1;
            matrix(i,j)=0;
            matrix(j,i)=0;
        end
    end
end
%reconnect the edge anticlockwisely for Node N
for i=N
    for j=1:m
        r=rand(1);
        if r<=p
            unconect=find(matrix(i,:)==0);
            r1=ceil(length(unconect)*rand(1));
            matrix(i,unconect(r1))=1;
            matrix(unconect(r1),i)=1;
            matrix(i,j)=0;matrix(j,i)=0;
        end
    end
end
for m=1:N
    % Avoid the self-loop
    matrix(m,m)=0;
end
A_adj=full(matrix);
% A_weight=add_weight(A_adj,N);
% plot_RanGraph(A_adj,0);



