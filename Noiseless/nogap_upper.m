% Numerical check gaps beween str = ii,iii, and iv
% verify no gap
clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;

valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
upper_res = zeros(3,10);
poolobj = parpool(10);
parfor itr = 1:10

    t = 0.3*itr;
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);

    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Ctheta = kron(Etheta,Etheta);
    % random wx
    rep = 1500;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end

    upper_res(:,itr) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2);
        SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
        SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];

end
delete(poolobj);
save nogap_upper_1500.mat upper_res;