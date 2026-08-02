function U = build_SchurWeyl_U(d)
    % BUILD_SCHURWEYL_U 构造 N=2, 任意局部维度 d 的 Schur-Weyl 基变换矩阵 U
    % U 的列向量由对称子空间基和反对称子空间基拼接而成
    
    dim_total = d^2;
    U = zeros(dim_total, dim_total); 
    
    col_idx = 1; 
 
    
    % a) 对角基 |ii>
    for i = 1:d
        row_idx = (i-1)*d + i;
        U(row_idx, col_idx) = 1;
        col_idx = col_idx + 1;
    end
    
    % b) 非对角对称基 1/sqrt(2) * (|ij> + |ji>)
    for i = 1:d
        for j = i+1:d
            row_idx_ij = (i-1)*d + j; % |ij> 的位置
            row_idx_ji = (j-1)*d + i; % |ji> 的位置
            
            U(row_idx_ij, col_idx) = 1/sqrt(2);
            U(row_idx_ji, col_idx) = 1/sqrt(2);
            col_idx = col_idx + 1;
        end
    end
    
    % c) 非对角反对称基 1/sqrt(2) * (|ij> - |ji>)
    for i = 1:d
        for j = i+1:d
            row_idx_ij = (i-1)*d + j; % |ij> 的位置
            row_idx_ji = (j-1)*d + i; % |ji> 的位置
            
            U(row_idx_ij, col_idx) = 1/sqrt(2);
            U(row_idx_ji, col_idx) = -1/sqrt(2);
            col_idx = col_idx + 1;
        end
    end
end