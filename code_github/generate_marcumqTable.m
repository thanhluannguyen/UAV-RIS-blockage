function generate_marcumqTable(K)

    fx = @(x) 1 - marcumq(sqrt(2*K), x);
    % fplot(fx, [0, 10])

    xx = linspace(0, 10, 1e4);
    yy = fx(xx);

    idx = min( find(yy > 1-10^(-5)) );

    xx = linspace(0, xx(idx), 1e2); 
    yy = fx(xx);

    marcumqTable = [yy.' xx.'];
    
    save('marcumqTable.mat', 'marcumqTable');
end