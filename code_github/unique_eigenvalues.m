function [sorted_unique_eigenvalues, multiplicities] = unique_eigenvalues(X)
    eigenvalues = eig(X);
    unique_eigenvalues = unique(eigenvalues);
    sorted_unique_eigenvalues = sort(unique_eigenvalues, 'descend');
    multiplicities = zeros(size(sorted_unique_eigenvalues));

    for i = 1:length(sorted_unique_eigenvalues)
        multiplicities(i) = sum(eigenvalues == sorted_unique_eigenvalues(i));
    end
end
