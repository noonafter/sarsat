function [max_f] = find_max_fft(data_freq)
%MAX_FFT_INTERP 利用fft内插求fft最大值对应的分数索引，0~Nfft-1
%   此处显示详细说明
Nfft = length(data_freq);
[midb1, max_idx] = max(abs(data_freq));
% interp
if max_idx == 1
   max_idx_left = Nfft;
   max_idx_right = 2;
elseif max_idx == Nfft
   max_idx_left = Nfft - 1;
   max_idx_right = 1;
else
   max_idx_left = max_idx - 1;
   max_idx_right = max_idx + 1;
end
midb0 = abs(data_freq(max_idx_left));
% midb1 = abs(data_freq(max_idx));
midb2 = abs(data_freq(max_idx_right));
midbdelta = 0.5 * (midb2 - midb0) / (2 * midb1 - midb2 - midb0);
max_f = mod(max_idx + midbdelta - 1,Nfft);

end

