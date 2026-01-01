clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
t = 0.1;
% true value, B = 1, theta = phi = pi/4;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

% state = RandomDensityMatrix(numVar);
% save state.mat state;
% choi matrix for the channel
% H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
% gamma = 0.1;
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
% res(:,idx) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1);
%     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2);
%     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
%     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];
% test res 
% res1 = SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1);
upper_res(:,idx) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2);
    SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
    SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];
end
delete(pool);
save upper_gamma66.mat upper_res;
% save upper1.mat res1;
% save upper2.mat res2;
% save upper3.mat res3;
% save upper4.mat res4;

% pool = parpool(40);
% parfor idx = 1:9
%     gamma = idx*0.1;
%     Kr1 = sqrt(eta)*eye(2); Kr2 = sqrt(1-eta)*[0 1;1 0];
%     Kr3 = [1 0; 0 sqrt(1-gamma)]; Kr4 = [0 sqrt(gamma);0 0];
% 
% 
%     U = expm(-1j*H*t);
%     % U = expm(-1j*theta1*Pauli(1)*t)*expm(-1j*theta2*Pauli(2)*t)*expm(-1j*theta3*Pauli(3)*t)*expm(-1j*cos(theta2)*RandomDensityMatrix(2)*t);
%     bigU = kron(eye(2),U);
%     bigKr1 = kron(eye(2),Kr1);
%     bigKr2 = kron(eye(2),Kr2);
%     bigKr3 = kron(eye(2),Kr3);
%     bigKr4 = kron(eye(2),Kr4);
%     Omega = [1;0;0;1]*[1 0 0 1];
%     Etheta = bigU*Omega*bigU';
%     Etheta = bigKr3*bigKr1*Etheta*bigKr1'*bigKr3'+bigKr3*bigKr2*Etheta*bigKr2'*bigKr3'+bigKr4*bigKr1*Etheta*bigKr1'*bigKr4'+bigKr4*bigKr2*Etheta*bigKr2'*bigKr4';
% %     disp(PartialTrace(Etheta,2,[2 2]))
% 
% 
% %     Etheta = Etheta + bigKr1*Etheta*bigKr1';
% %     Etheta = Etheta + bigKr2*Etheta*bigKr2';
% %     Etheta = Etheta + bigKr3*Etheta*bigKr3';
% %     Etheta = Etheta + bigKr4*Etheta*bigKr4';
%     Ctheta = kron(Etheta,Etheta);
% 
%     % random wx
%     rep = 100;
%     wx = zeros(numVar+1,rep);
%     for i = 1:rep
%         wx(:,i) = RandomStateVector(numVar+1,1);
%     end
%     % res(:,idx) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1);
%     %     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2);
%     %     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 3);
%     %     SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 4)];
%     % test res 
%     res(:,idx) = [SDP_upper(N, wx, Ctheta, theta1, valth1, 2);
%          SDP_upper(N, wx, Ctheta, theta1, valth1, 3)];
% end
% delete(pool);
% save res_100.mat res;