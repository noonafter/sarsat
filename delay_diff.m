function [data_out] = delay_diff(data_in, delay)
%DELAY_DIFF 此处显示有关此函数的摘要
%   此处显示详细说明
data_out = data_in .* conj([zeros(1,delay) data_in(1:end-delay)]);

end

