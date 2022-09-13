function [mseq_i,mseq_q] = mseq_gene(mseq_len)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% mseq_len = 38400;
if(nargin < 1)
    mseq_len = 38400+64;
end
mseq_i = zeros(1,mseq_len);
mseq_q = zeros(1,mseq_len);

g = zeros(1,23);
g(5) = 1;
g(23) = 1;

reg_i = [zeros(1,22) 1]; % initial state
reg_q = [0 0 1 ...
    1 0 1 0 ...
    1 1 0 0 ...
    0 0 0 1 ...
    1 1 1 1 ...
    1 1 0 0]; % initial state

mseq_i(1) = reg_i(end);
mseq_q(1) = reg_q(end);
for index = 2:mseq_len
tmp_i = mod(g * reg_i',2);
tmp_q = mod(g * reg_q',2);
reg_i = [tmp_i reg_i(1:22)];
reg_q = [tmp_q reg_q(1:22)];
mseq_i(index) = reg_i(23);
mseq_q(index) = reg_q(23);
end


end

