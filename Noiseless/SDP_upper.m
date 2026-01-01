function [res,X,newXk] = SDP_upper(N, wx, Ctheta, var, varvals, str)
% str stands for types of strategies. 1 for i, 2 for ii, 3 for iii and 4
% for iv
% wx is a matrix whose ith column is the vector |w_x>
    numvec = size(wx,2);
    numVar = size(var,2); % require var to be a row vector
    dimCalC = numVar + 1;
    dimO = 2^N;
    dimIO = 4^N;
	tildeW = blkdiag(0,eye(3));
    newXk = zeros(dimIO,dimIO,numvec);
%     cvx_solver mosek;
    cvx_begin sdp 
        variable Xk(dimIO,dimIO,numvec) complex semidefinite;
        variable X1(dimIO,dimIO) complex semidefinite;
        variable X2(dimIO,dimIO) complex semidefinite;
        expression X(dimIO,dimIO);
        expression judDelta(numVar,numVar);
        expression bigX(dimCalC*dimIO, dimCalC*dimIO);
        bigX = 0;
        
        for idx = 1:numvec
            bigX = bigX + kron(wx(:,idx)*wx(:,idx)',Xk(:,:,idx));
        end
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
                    )*bigX);
            end
        end
%         X = zeros(dimIO,dimIO);
        X = 0;
        for idx = 1:numvec
            X = X + wx(1,idx)^2*Xk(:,:,idx);
        end
%       restrictions
        judDelta == eye(numVar);
        
        trace(X) == dimO;
        if str == 1
%             parallel
            myoperation(X,2*(1:N),repmat(2, 1, 2*N)) == X;
        elseif str == 2
%             sequential only implement the N = 2 case now
            myoperation(X,4,[2 2 2 2]) - myoperation(X,[3 4],[2 2 2 2]) + myoperation(X,[2 3 4],[2 2 2 2]) == X;
        elseif str == 3
            X == X1 + X2;
            myoperation(X1,4,[2 2 2 2]) - myoperation(X1,[3 4],[2 2 2 2]) + myoperation(X1,[2 3 4],[2 2 2 2]) == X1;
            myoperation(X2,2,[2 2 2 2]) - myoperation(X2,[1 2],[2 2 2 2]) + myoperation(X2,[1 2 4],[2 2 2 2]) == X2;
        elseif str == 4
            myoperation(X,[1 2 4],[2 2 2 2]) - myoperation(X,[1 2],[2 2 2 2]) ...
                + myoperation(X,[2 3 4],[2 2 2 2]) - myoperation(X,[3 4],[2 2 2 2]) ...
                + myoperation(X,2,[2 2 2 2]) + myoperation(X,4,[2 2 2 2]) - myoperation(X,[2 4],[2 2 2 2]) == X;
        end
%       objective function
        minimize(real(trace(kron(tildeW,double(subs(Ctheta,var,varvals)))*bigX)))
    cvx_end
    res = cvx_optval;
    for i = 1:numvec
        newXk(:,:,i) = wx(1,i)^2*Xk(:,:,i);
    end
end