clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
t = 0.1;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
pool = parpool(10);
upper_res = zeros(3,9);
parfor idx = 1:9

    gamma = 0.66 + idx*0.01;
Kr3 = [1 0; 0 sqrt(1-gamma)]; Kr4 = [0 sqrt(gamma);0 0];
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
bigKr3 = kron(eye(2),Kr3);
bigKr4 = kron(eye(2),Kr4);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Etheta = bigKr3*Etheta*bigKr3' + bigKr4*Etheta*bigKr4';
Ctheta = kron(Etheta,Etheta);
% random wx
rep = 1500;
wx = zeros(numVar+1,rep);
for i = 1:rep
    wx(:,i) = RandomStateVector(numVar+1,1);
end
upper_res(:,idx) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2);
    SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
    SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];
end
delete(pool);
save upper_gamma66.mat upper_res;