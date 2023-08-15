function out = getTau(K, L, blockedProb)
    marcumqTable = cell2mat(struct2cell(load('marcumqTable.mat')));
    
    [~, ind] = min( abs(blockedProb-marcumqTable(:, 1)) );
    
    out = marcumqTable(ind, 2)^2 *(L/(2*(K+1)));
end