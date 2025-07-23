% compute the trace of a product of a matrix x 
% with a tensor product of Toeplitz matrices. 

function TrX = TrFind(x, n, k)

    %producing the required tensor product of Toeplitz matrices
    TrX = 0;
    y = 1;
    for i = 1 : length(k)
        if(abs(k(i)) > n(i)) % (offset from the diagonal)
            y = kron(y, zeros(n(i)+1, n(i)+1));
        else
            y = kron(y,diag(ones(n(i)+1 - abs(k(i)), 1), k(i))); % toeplitz matrix is created
        end
    end
    TrX = TrX + trace(y' * x); % outputs the value of Tr_k^n(X)

end
