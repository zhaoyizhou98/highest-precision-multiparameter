clear all; clc;

N = 2;
syms theta1; syms theta2; syms theta3;
numVar = 3;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
lower_res = zeros(1,10);


t = 0.05*1+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,1) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*2+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,2) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*3+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,3) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*4+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,4) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*5+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,5) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*6+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,6) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*7+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,7) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*8+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,8) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*9+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,9) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

t = 0.05*10+2.5;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);

save lower_t2dot5to3.mat lower_res;