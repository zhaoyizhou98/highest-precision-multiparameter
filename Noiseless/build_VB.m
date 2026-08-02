function V_B = construct_VB(d_B, n)
    % 计算对称子空间的维度 s_n
    s_n = nchoosek(d_B + n - 1, n);
    
    % 如果 n=0，直接返回 1
    if n == 0
        V_B = 1;
        return;
    end
    
    % 生成所有非递减的多重指标（多重集）
    M = generate_multiset(d_B, n);
    total_dim = d_B^n;
    V_B = zeros(total_dim, s_n);
    base = d_B.^(0:n-1); % 用于计算线性索引的基数
    
    % 遍历每个多重指标
    for idx = 1:s_n
        I = M(idx, :);
        % 生成所有不同的排列
        P = unique(perms(I), 'rows');
        K = size(P, 1);
        v = zeros(total_dim, 1);
        
        % 对每个排列计算其在张量积空间中的索引
        for j = 1:K
            p = P(j, :);
            ind = 1 + sum((p - 1) .* base);
            v(ind) = v(ind) + 1;
        end
        v = v / sqrt(K); % 归一化
        V_B(:, idx) = v;
    end
end

function M = generate_multiset(d, n)
    M = generate_multiset_from(1, d, n);
end

function M = generate_multiset_from(start, d, n)
    if n == 0
        M = [];
    elseif n == 1
        M = (start:d)';
    else
        M = [];
        for i = start:d
            M_sub = generate_multiset_from(i, d, n-1);
            M = [M; [i * ones(size(M_sub, 1), 1) M_sub]];
        end
    end
end