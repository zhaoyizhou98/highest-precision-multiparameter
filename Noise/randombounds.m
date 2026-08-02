% potential universality across different noise models
% random channels

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
data = load("RandomChannels4Ele.mat");
K = data.K;
% % % % % % % % % % % % % % % % % % % % % % % 
% lower bounds
% % % % % % % % % % % % % % % % % % % % % % % 
for col = 1:8
    Kr1 = K{col}(:,:,1); Kr2 = K{col}(:,:,2); Kr3 = K{col}(:,:,3); Kr4 = K{col}(:,:,4);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr1 = kron(eye(2),Kr1);
    bigKr2 = kron(eye(2),Kr2);
    bigKr3 = kron(eye(2),Kr3);
    bigKr4 = kron(eye(2),Kr4);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2' + bigKr3*Etheta*bigKr3' + bigKr4*Etheta*bigKr4';
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

    lower_res = [SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},1);
        SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},2);
        SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},3);
        SDP_lower_Num(2,N,Ctheta,{Ctheta1,Ctheta2,Ctheta3},4)];
        
     save(sprintf('4Ele_random_lower_col_%d.mat', col), 'lower_res');
        % save random_lower.mat lower_res;
end

%%
% random wx
rep = 3500;
wx = zeros(numVar+1,rep);
for i = 1:rep
    wx(:,i) = RandomStateVector(numVar+1,1);
end
data = load("RandomChannels4Ele.mat");
K = data.K;

% % % % % % % % % % % % % % % % % % % % % % % 
% upper bounds
% % % % % % % % % % % % % % % % % % % % % % % 
for idx = 3:8
    Kr1 = K{idx}(:,:,1); Kr2 = K{idx}(:,:,2); Kr3 = K{idx}(:,:,3); Kr4 = K{idx}(:,:,4);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    bigKr1 = kron(eye(2),Kr1);
    bigKr2 = kron(eye(2),Kr2);
    bigKr3 = kron(eye(2),Kr3);
    bigKr4 = kron(eye(2),Kr4);
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Etheta = bigKr1*Etheta*bigKr1' + bigKr2*Etheta*bigKr2' + bigKr3*Etheta*bigKr3' + bigKr4*Etheta*bigKr4';
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

    upper_res = zeros(3,1);
    pool = parpool(3);
    parfor str = 1:3
        upper_res(str) = SDP_upper_Num(N, wx, Ctheta, {Ctheta1,Ctheta2,Ctheta3}, str+1);
    end
    delete(pool);
    save(sprintf('4Ele_random_upper_idx_%d.mat', idx), 'upper_res');
end
%%
clear all; clc;
for idx = 1:8
    data1 = load(sprintf('4Ele_random_upper_idx_%d.mat', idx));
    data2 = load(sprintf('4Ele_random_lower_col_%d.mat', idx));
    data2.lower_res(1:3,1) - data1.upper_res(:,1)
end
%%
clear all; clc;
for idx = 1:8
    K{idx} = myrandomChannel(2,4);
end

save RandomChannels4Ele.mat K;