function [A,N] = read_large_sparse_network(path)

    fd = fopen(path,'rt');
    formatSpec = '%d %d';
    sizeA = [2 Inf];

    Y = fscanf(fd, formatSpec, sizeA);
    Y = Y';
    fclose(fd);
    N = max(max(Y(:,1)),max(Y(:,2)));
    t = size(Y, 1);
    o = ones(t,1);
    A = sparse(Y(:,1),Y(:,2),o,N,N);
    
    clear Y o;
end

