function N = id_neigh(Adj)
% identify set of neighbors (for busses and agents)
% W. Ananduta
% 30/09/2020

for i = 1:size(Adj,1)
    [~,N{i}] = find(Adj(i,:));
end

end