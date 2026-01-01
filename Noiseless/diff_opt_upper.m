% t = 2.5 to 3
clear all; clc;
% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
% true value, B = 1, theta = phi = pi/4;
valth1 = 1/2; valth2 = 1/2; valth3 = sqrt(2)/2;

H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);

upper_res = zeros(2,10);
poolobj = parpool(10);
parfor itr = 1:10
    t = 0.05*itr+2.5;
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Ctheta = kron(Etheta,Etheta);
% random wx
    rep = 125;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end

    % another method based on works of Hayashi and Ouyang, which is not random
    % repn = 5;
    % wx = generateDn_d(repn,numVar)';
    
    upper_res(:,itr) = [SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1);
        SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 2)];
end
delete(poolobj);
save 5run_t2dot5to3.mat upper_res;
%%
% compare with analytical bounds (Application 1)
clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
% true value, B = 1, theta = phi = pi/4;

len = 20;
data = zeros(len+1,len+1);

poolobj = parpool(20);
t = 3;
temp_data = zeros(1,(1+len)^2);
parfor itr = 1:(1+len)^2
    disp("Current loop: "+itr);
    xaxis = linspace(0,1,len+1); yaxis = linspace(1,0,len+1);
    [X, Y] = meshgrid(xaxis,yaxis);
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    valth1 = X(row,col); valth2 = Y(row,col);
    if valth1^2 + valth2^2 > 1
        continue;
    end
    
    valth3 = sqrt(1-valth1^2-valth2^2);

    H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Ctheta = kron(Etheta,Etheta);
    % random wx
    rep = 700;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end
    % another method based on works of Hayashi and Ouyang, which is not random
    % repn = 9;
    % wx = generateDn_d(repn,numVar)';
    temp_data(itr) = SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1); 
end

for itr = 1:(1+len)^2
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    data(row, col) = temp_data(itr);
end
save colormap_data_20.mat data;

%%
clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
% true value, B = 1, theta = phi = pi/4;

len = 20;
data = zeros(len+1,len+1);

poolobj = parpool(20);
t = 3;
temp_data = zeros(1,(1+len)^2);
parfor itr = 1:(1+len)^2
    disp("Current loop: "+itr);
    xaxis = linspace(0,sqrt(2),len+1); yaxis = linspace(sqrt(2),0,len+1);
    [X, Y] = meshgrid(xaxis,yaxis);
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    valth1 = X(row,col); valth2 = Y(row,col);
    if valth1^2 + valth2^2 > 2
        continue;
    end
    
    valth3 = sqrt(2-valth1^2-valth2^2);

    H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Ctheta = kron(Etheta,Etheta);
    % random wx
    rep = 700;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end
    % repn = 9;
    % wx = generateDn_d(repn,numVar)';
    temp_data(itr) = SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1); 
end

for itr = 1:(1+len)^2
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    data(row, col) = temp_data(itr);
end
save colormap_data_20_norm2.mat data;

%%
clear all; clc;

% N = 2 case, i.e., two signal channels
N = 2;
% init
syms theta1; syms theta2; syms theta3;
numVar = 3;
% true value, B = 1, theta = phi = pi/4;

len = 20;
data = zeros(len+1,len+1);

poolobj = parpool(20);
t = 3;
temp_data = zeros(1,(1+len)^2);
parfor itr = 1:(1+len)^2
    disp("Current loop: "+itr);
    xaxis = linspace(0,sqrt(3),len+1); yaxis = linspace(sqrt(3),0,len+1);
    [X, Y] = meshgrid(xaxis,yaxis);
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    valth1 = X(row,col); valth2 = Y(row,col);
    if valth1^2 + valth2^2 > 3
        continue;
    end
    
    valth3 = sqrt(3-valth1^2-valth2^2);

    H = theta1*Pauli(1) + theta2*Pauli(2) + theta3*Pauli(3);
    U = expm(-1j*H*t);
    bigU = kron(eye(2),U);
    
    Omega = [1;0;0;1]*[1 0 0 1];
    Etheta = bigU*Omega*bigU';
    Ctheta = kron(Etheta,Etheta);
    % random wx
    rep = 700;
    wx = zeros(numVar+1,rep);
    for i = 1:rep
        wx(:,i) = RandomStateVector(numVar+1,1);
    end
    % repn = 9;
    % wx = generateDn_d(repn,numVar)';
    temp_data(itr) = SDP_upper(N, wx, Ctheta, [theta1,theta2,theta3], [valth1 valth2 valth3], 1); 
end

for itr = 1:(1+len)^2
    row = ceil(itr/21); col = mod(itr-1,21)+1;
    data(row, col) = temp_data(itr);
end
save colormap_data_20_norm3.mat data;
%%
function [row, col] = find_pos_by_dimension(num, N)

    if num <= 0 || N <= 0
        error('Input number and dimension of the matrix must be positive integers');
    end
    
    max_val = N * (N + 1) / 2;
    
    if num > max_val
        error('Number %d exceeds %d dimension (max %d)', num, N, max_val);
    end

    current_start_val = 1; 
    
    for c = 1:N
        col_len = N - c + 1;
        
        col_end_val = current_start_val + col_len - 1;
        
        if num >= current_start_val && num <= col_end_val
            col = c;
            
            offset = num - current_start_val;
            row = c + offset;
            
            return;
        end
        
        current_start_val = col_end_val + 1;
    end
end