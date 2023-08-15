function M = generate_partitions(N)
    M = cell(0);
    memo = cell(N, N);
    for L = 1:N
        partial_partition = cell(1, L);
        partitions = generate_partitions_recursive(N, L, 1, partial_partition, memo);
        M = [M; partitions];
    end
end

function partitions = generate_partitions_recursive(N, L, min_value, partial_partition, memo)
    if L == 1
        if N >= min_value
            partial_partition{1} = N;
            partitions = {cell2mat(partial_partition)};
        else
            partitions = {};
        end
    else
        if ~isempty(memo{N, L}) && min_value == 1
            partitions = memo{N, L};
        else
            partitions = {};
            for i = min_value:N-1
                partial_partition{L} = i;
                sub_partitions = generate_partitions_recursive(N-i, L-1, i, partial_partition, memo);
                if ~isempty(sub_partitions)
                    partitions = [partitions; sub_partitions];
                end
            end
            if min_value == 1
                memo{N, L} = partitions;
            end
        end
    end
end
