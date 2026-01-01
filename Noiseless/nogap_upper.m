clear all; clc;

N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
t = 3;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);

U = expm(-1j*H*t);
bigU = kron(eye(2),U);

Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
% random wx
rep = 1000;
wx = zeros(numVar+1,rep);
for i = 1:rep
    wx(:,i) = RandomStateVector(numVar+1,1);
end

upper_res = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
    SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];

save nogap_upper.mat upper_res;