clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
t = 0.1;
% true value, B = 1, theta = phi = pi/4;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);

value = zeros(3,9);

for col = 1:9
    gamma = 0.66+col*0.01;
    Kr3 = [1 0; 0 sqrt(1-gamma)]; Kr4 = [0 sqrt(gamma);0 0];
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr3 = kron(eye(2),Kr3);
    bigKr4 = kron(eye(2),Kr4);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr3*Etheta*bigKr3' + bigKr4*Etheta*bigKr4';
    Ctheta = kron(Etheta,Etheta);
    
    pool = parpool(4);
    parfor idx = 1:3
    value(idx,col) = SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],idx);
    end
    delete(pool);
end

save lower_gamma66.mat value;