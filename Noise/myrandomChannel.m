% function K = myrandomChannel(dim)
% % Randomly generate quantum channel by randomly generating Unitary matrices
% % Input: dim stands for the input and output dimension of the quantum
% % channel
% % Output: K contains Kraus operators. Get one Kraus operators by K(:,:,i)
% 
% K = zeros(dim,dim,dim);
% U = RandomUnitary(dim*dim);
% for i = 1:dim
%     K(:,:,i) = U((i-1)*dim+1:i*dim,1:dim);
% end
% end

function K = myrandomChannel(dim,numKraus)
% Random quantum channel with numKraus Kraus operators
%
% dim       : input/output dimension
% numKraus  : number of Kraus operators
%
% K(:,:,i) is the i-th Kraus operator

K = zeros(dim,dim,numKraus);

U = RandomUnitary(dim*numKraus);

for i = 1:numKraus
    K(:,:,i) = U((i-1)*dim+1:i*dim,1:dim);
end

end