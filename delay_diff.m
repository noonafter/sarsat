function [data_out] = delay_diff(data_in, delay)
%DELAY_DIFF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
data_out = data_in .* conj([zeros(1,delay) data_in(1:end-delay)]);

end

