function res = SDP_lower(n,N, Ctheta, var, varvals, str)
% str stands for types of strategies. 1 for i, 2 for ii, 3 for iii and 4
% for iv
    numVar = size(var,2); % require var to be a row vector
    dimCalC = numVar + 1;
    dimO = 2^N;
    dimIO = 4^N;
    tildeW = blkdiag(0,eye(3));
    d = 4;
    S = zeros(d^2, d^2);
    for i = 1:d
        for j = 1:d
            eij = zeros(d,d); eij(i,j) = 1;       
            eji = zeros(d,d); eji(j,i) = 1;       
            S = S + kron(eij, eji);
        end
    end
%     cvx_solver mosek
    cvx_begin sdp quiet
        variable Ykn(dimCalC^2*dimIO,dimCalC^2*dimIO) complex semidefinite;
        variable PartialYk1(dimIO,dimIO) complex semidefinite;
        variable PartialYk2(dimIO,dimIO) complex semidefinite;
        expression Yk(dimCalC*dimIO,dimCalC*dimIO);
        Yk = PartialTrace(Ykn,1:n-1,[repmat(dimCalC, 1, n), dimIO]);
        expression judDelta(numVar,numVar);
        for i = 1:numVar
            for j = 1:numVar
                tmpmat = zeros(dimCalC,dimCalC);
                tmpmat(i+1,1) = 1/2;
                tmpmat(1,i+1) = 1/2;
                judDelta(i,j) = trace(kron(tmpmat, ... 
                    double( ... 
                    subs( ... 
                    diff(Ctheta,var(j)),var,varvals) ... 
                    ) ...
                    )*Yk);
            end
        end
        Ykn == kron(S,eye(dimIO))*Ykn*kron(S,eye(dimIO))';
        expression PartialYk(dimIO,dimIO);
        tmpmat = zeros(dimCalC,dimCalC);
        tmpmat(1,1) = 1;
        PartialYk = PartialTrace(kron(tmpmat,eye(dimIO))*Yk,[1],[dimCalC,repmat(2, 1, 2*N)]);
% % % % % % % %       restrictions
        judDelta == eye(numVar);
%         second restrictions, very difficult
        trace(PartialYk) == dimO;
        if str == 1
%             parallel strategies
            myoperation(PartialYk,2*(1:N),repmat(2, 1, 2*N)) == PartialYk;
        elseif str == 2
            myoperation(PartialYk,4,[2 2 2 2]) - myoperation(PartialYk,[3 4],[2 2 2 2]) + myoperation(PartialYk,[2 3 4],[2 2 2 2]) == PartialYk;
        elseif str == 3
            PartialYk == PartialYk1 + PartialYk2;
            myoperation(PartialYk1,4,[2 2 2 2]) - myoperation(PartialYk1,[3 4],[2 2 2 2]) + myoperation(PartialYk1,[2 3 4],[2 2 2 2]) == PartialYk1;
            myoperation(PartialYk2,2,[2 2 2 2]) - myoperation(PartialYk2,[1 2],[2 2 2 2]) + myoperation(PartialYk2,[1 2 4],[2 2 2 2]) == PartialYk2;
        elseif str == 4
            myoperation(PartialYk,[1 2 4],[2 2 2 2]) - myoperation(PartialYk,[1 2],[2 2 2 2]) ...
                + myoperation(PartialYk,[2 3 4],[2 2 2 2]) - myoperation(PartialYk,[3 4],[2 2 2 2]) ...
                + myoperation(PartialYk,2,[2 2 2 2]) + myoperation(PartialYk,4,[2 2 2 2]) - myoperation(PartialYk,[2 4],[2 2 2 2]) == PartialYk;
        end
        for idx = 1:n
            PartialTranspose(Ykn,1:idx,[repmat(dimCalC, 1, n), dimIO]) >= 0;
        end
% % % % % % % %     objective function 
        minimize(real(trace(kron(tildeW,double(subs(Ctheta,var,varvals)))*Yk)))
    cvx_end
    res = cvx_optval;
end