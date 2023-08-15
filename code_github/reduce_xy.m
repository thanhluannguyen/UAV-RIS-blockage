function [x_new, y_new] = reduce_xy(x, y, n, mode)
    if strcmp(mode, 'Log')
        max_idx = floor( log(length(x))/log(n) );
        idx = n.^(1:max_idx);
    elseif strcmp(mode, 'Linear')
        idx = 1:n:length(x);
    end

    x_new = x(idx);
    y_new = y(idx);
end