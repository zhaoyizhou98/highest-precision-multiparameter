% SDP lower bound with symmetry
% N = 2 case
% can only be applied to str = i and iv
function res = Sym_SDP_lower_Num(n,N, Ctheta, var, str)
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
    sn = nchoosek(numVar+n,n);
    % input and output are qubits
    d11 = d*(d-1)/2;
    d20 = d*(d+1)/2;
    V = build_VB(dimCalC,n);
    U = build_SchurWeyl_U(d);
%   cvx_solver mosek
    cvx_begin sdp
        % % % % % % % %         variables  
        variable PartialYk1(dimIO,dimIO) complex semidefinite;
        variable PartialYk2(dimIO,dimIO) complex semidefinite;
        variable Z11(sn*d11,sn*d11) complex semidefinite;
        variable Z20(sn*d20,sn*d20) complex semidefinite;
        % expression Yk(dimCalC*dimIO,dimCalC*dimIO);
        % expression Ykn(dimCalC^2*dimIO,dimCalC^2*dimIO);
        expression Yk(dimCalC*dimIO,dimCalC*dimIO);
        expression Ykn(dimCalC^2*dimIO,dimCalC^2*dimIO);
        UV_full = kron(U,V);
        Zbig = [Z20, zeros(size(Z20,1), size(Z11,2));
        zeros(size(Z11,1), size(Z20,2)), Z11];
        Ykn = UV_full*Zbig*UV_full';
        Yk = Swap(PartialTrace(Ykn,3:n+1,[dimIO,repmat(dimCalC, 1, n)]),[1,2],[dimIO,dimCalC]);
        for k = 1:n
            PartialTranspose(Ykn,[2:k+1],[dimIO,repmat(dimCalC, 1, n)]) >= 0;
        end
        expression judDelta(numVar,numVar);
        for i = 1:numVar
            for j = 1:numVar
                tmpmat = zeros(dimCalC,dimCalC);
                tmpmat(i+1,1) = 1/2;
                tmpmat(1,i+1) = 1/2;
                judDelta(i,j) = trace(kron(tmpmat,var{j})*Yk);
            end
        end
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
        elseif str == 4
            myoperation(PartialYk,[1 2 4],[2 2 2 2]) - myoperation(PartialYk,[1 2],[2 2 2 2]) ...
                + myoperation(PartialYk,[2 3 4],[2 2 2 2]) - myoperation(PartialYk,[3 4],[2 2 2 2]) ...
                + myoperation(PartialYk,2,[2 2 2 2]) + myoperation(PartialYk,4,[2 2 2 2]) - myoperation(PartialYk,[2 4],[2 2 2 2]) == PartialYk;
        end
% % % % % % % %     objective function 
        minimize(real(trace(kron(tildeW,Ctheta)*Yk)))
    cvx_end
    res = cvx_optval;
end