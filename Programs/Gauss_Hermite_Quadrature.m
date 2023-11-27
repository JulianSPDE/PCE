function [x, w] = Gauss_Hermite_Quadrature(n)
    % n: number of nodes and weights
    
    % Construct the companion matrix
    C = zeros(n);
    C(2:n+1:n*(n-1)) = sqrt(1:n-1)/sqrt(2);
    C = C + C';
    
    % Compute the eigenvalues and eigenvectors of the companion matrix
    [V, D] = eig(C);
    
    % Sort the eigenvalues and eigenvectors in ascending order
    [lambda, indices] = sort(diag(D));
    V = V(:, indices);
    
    % Compute the nodes and weights
    x = lambda;
    w = sqrt(pi) * V(1, :).^2;
end