% check more standard noise models
clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
var = [theta1,theta2,theta3];
numVar = 3;
t = 1;
% true value, B = 1, theta = phi = pi/4;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;
varvals = [valth1,valth2,valth3];

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);


%%
% % % % % % % % % % % % % % % % % % % % % % % 
% dephasing noise
% % % % % % % % % % % % % % % % % % % % % % % 
for cnt = 1:9
    p = cnt*0.1;
    Kr1 = sqrt(p)*eye(2); Kr2 = sqrt(1-p)*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr1 = kron(eye(2),Kr1);
    bigKr2 = kron(eye(2),Kr2);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2';
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


% % 
% lower bound
% % 

    lower_res = zeros(4,1);
    for idx = 1:4
        lower_res(idx,1) = SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},idx);
    end
    
    save(sprintf('dephasing_lower_%d.mat', cnt), 'lower_res');
end
%%
% % 
% upper bound
% % 
p = 0.1;
Kr1 = sqrt(p)*eye(2); Kr2 = sqrt(1-p)*Pauli(3);
U = expm(-1j*H*t);
bigU = kron(eye(2),U);
bigKr1 = kron(eye(2),Kr1);
bigKr2 = kron(eye(2),Kr2);
Omega = [1;0;0;1]*[1 0 0 1];
Etheta = bigU*Omega*bigU';
Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2';
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

% upper_res = zeros(3,1);
rep = 6000;
wx = zeros(numVar+1,rep);
for i = 1:rep
    wx(:,i) = RandomStateVector(numVar+1,1);
end
% upper_res = [SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 2);
%     SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 3);
%     SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 4)];
pool = parpool(2);
upper_res = zeros(2,1);
parfor itr = 1:2
    upper_res(itr,1) = SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, itr+2);
end
delete(pool);
% upper_res = [SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, 3);
%     SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, 4)];
% save(sprintf('dephasing_upper_rep4500_p_%d.mat', cnt), 'upper_res');
save dephasing_upper_p1_rep6000.mat upper_res;
%%
% % % % % % % % % % % % % % % % % % % % % % % 
% depolarizing noise
% % % % % % % % % % % % % % % % % % % % % % % 
for cnt = 1:13
    p = cnt*0.1;
    Kr1 = sqrt(1-3*p/4)*eye(2); Kr2 = sqrt(p)/2*Pauli(1); Kr3 = sqrt(p)/2*Pauli(2); Kr4 = sqrt(p)/2*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr1 = kron(eye(2),Kr1);
    bigKr2 = kron(eye(2),Kr2);
    bigKr3 = kron(eye(2),Kr3);
    bigKr4 = kron(eye(2),Kr4);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2'+bigKr3*Etheta*bigKr3'+bigKr4*Etheta*bigKr4';
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
    
    
    % % 
    % lower bound
    % % 
    
    lower_res = zeros(4,1);
    for idx = 1:4
        lower_res(idx,1) = SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},idx);
    end
    save(sprintf('depolarizing_lower_%d.mat', cnt), 'lower_res');
end
%%
% % 
% upper bound
% % 
for cnt  = 1:13
    p = cnt*0.1;
    Kr1 = sqrt(1-3*p/4)*eye(2); Kr2 = sqrt(p)/2*Pauli(1); Kr3 = sqrt(p)/2*Pauli(2); Kr4 = sqrt(p)/2*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr1 = kron(eye(2),Kr1);
    bigKr2 = kron(eye(2),Kr2);
    bigKr3 = kron(eye(2),Kr3);
    bigKr4 = kron(eye(2),Kr4);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2'+bigKr3*Etheta*bigKr3'+bigKr4*Etheta*bigKr4';
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
    % upper_res = zeros(3,1);
    rep = 4500;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end
    pool = parpool(2);
    upper_res = zeros(2,1);
    parfor itr = 1:2
        upper_res(itr,1) = SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, itr+2);
    end
    delete(pool);
% upper_res = [SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, 3);
%     SDP_upper_Num(N, wx, Ctheta,  {Ctheta1,Ctheta2,Ctheta3}, 4)];
% save(sprintf('dephasing_upper_rep4500_p_%d.mat', cnt), 'upper_res');
    
    % upper_res = [SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 2);
    %     SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 3);
    %     SDP_upper_Num(N, wx, Ctheta, [theta1,theta2,theta3], {Ctheta1,Ctheta2,Ctheta3}, 4)];
    save(sprintf('depolarizing_upper_rep4500_p_%d.mat', cnt), 'upper_res');
end