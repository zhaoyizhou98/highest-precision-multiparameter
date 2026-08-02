clear all; clc;

N = 2;
syms theta1; syms theta2; syms theta3;
var = [theta1,theta2,theta3];
numVar = 3;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;
varvals = [valth1,valth2,valth3];

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
lower_res = zeros(1,14);
sym_lower_res = 0;

t = 0.01*1+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
Ctheta1 = double( ... 
                subs( ... 
                diff(Ctheta,var(1)),var,varvals) ... 
                );
Ctheta2 = double( ... 
                subs( ... 
                diff(Ctheta,var(2)),var,varvals) ... 
                );
Ctheta3 = double( ... 
                subs( ... 
                diff(Ctheta,var(3)),var,varvals) ... 
                );
Ctheta = double(subs(Ctheta,var,varvals));

% lower_res(1,1) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
sym_lower_res = Sym_SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},1);
save sym_lower.mat sym_lower_res;

%%
t = 0.01*2+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,2) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("2")

t = 0.01*3+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,3) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("3")

t = 0.01*4+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,4) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("4")

t = 0.01*5+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,5) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("5")

t = 0.01*6+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,6) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("6")

t = 0.01*7+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,7) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("7")

t = 0.01*8+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,8) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("8")

t = 0.01*9+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,9) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("9")

t = 0.01*10+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("10")

t = 0.01*11+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("11")

t = 0.01*12+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("12")

t = 0.01*13+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("13")

t = 0.01*14+2.99;
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Ctheta = kron(Etheta,Etheta);
lower_res(1,10) =  SDP_lower(2,N,Ctheta,[theta1,theta2,theta3],[valth1 valth2 valth3],1);
disp("14")

save lower_t3to3dot13.mat lower_res;