function mu = mutual_coherence(A)
% First l2-normalize the columns of A. 
% Compute the dot products using A*A (Gram matrix)
A_S=bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));
mu=max(max(triu(transpose(A_S)*A_S,1)));
end

