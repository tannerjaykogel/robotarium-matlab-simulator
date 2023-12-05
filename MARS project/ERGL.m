function [L] = ERGL(n,p)
% creates a graph laplacian for an ER random graph with n nodes and
% probability p

L = zeros(n);   % initialize laplacian

% get edges that exist
for i = 1:n
    for j = 1:n
        if rand <= p && i ~= j
            % input edge to laplacian
            L(i,j) = -1; L(j,i) = -1;
        end
    end
end

% input degree of each node on diagonal
L = L - diag(sum(L,2));

end