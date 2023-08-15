clear all;

max_no_record = 40;
record = cell(max_no_record, 1, 1);
for irecord = 1:max_no_record
    M = generate_partitions(irecord);
    
    coeff = 0;
    rec_coeff = [];
    rec_eigs = [];
    rec_ms = [];
    for iM = 1:length(M)
        [eigs, ms] = unique_eigenvalues(diag(M{iM}));

        coeff = prod(1./factorial(ms)) / prod(factorial(M{iM}));
        
        n_zeros = (irecord-1)-length(eigs);
        eigs = [eigs; ones(n_zeros, 1)];

        n_zeros = (irecord-1)-length(ms);
        ms = [ms; zeros(n_zeros, 1)];
        
        rec_coeff = [rec_coeff coeff];
        rec_eigs = [rec_eigs eigs];
        rec_ms = [rec_ms ms];
    end
    record{irecord, 1} = rec_coeff;
    record{irecord, 2} = rec_eigs;
    record{irecord, 3} = rec_ms;
end

save('record.mat', 'record');