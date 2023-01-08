function [A,edges,nset] = plotbusV(bus,Lnbr)

edges = Lnbr(:,2:3);
nb = size(bus,1);
s = edges(:,1); t = edges(:,2);
G = digraph(s,t);

% find all terminal nodes
I = full(incidence(G)); 
target=[];
for i = 1 : nb
if sum(I(i,:))>0
    target = [target;i];
end
end

[TR,D] = shortestpathtree(G,1,target,'OutputForm','cell');

nset = size(TR,1);
A={};
for i = 1 : nset
    k=1;
    for j = 1  : size(TR{i},2)
        id = find(TR{i}(1,j,1)==edges(:,1));
        n_id= size(id,1);
        for J = 1 : n_id
            if  ~isempty(find(edges(id(J,1),2)==TR{i}(1,:,1)))
                A{i}(:,k,:)=id(J,1);
                k=k+1;
            end
        end
    end  
end


for i = 2 : nset
    A{i} = setdiff(A{i},A{1});  
end
 

