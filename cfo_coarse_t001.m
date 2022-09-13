function [data_out, sig_pos, CFO_est] = cfo_coarse_t001(data_in, fsa, Nfft, thd)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
pass_thd_count = 0;
cfo_dec = 0;
sig_pos = -1;
% tmp = [];
% CFO estimation
loop_num = floor(length(data_in) / Nfft);
for idx = 1:loop_num
    data_cnt = data_in((idx-1)*Nfft + [1:Nfft]);
    data_freq =  fft(data_cnt, Nfft);
    max_f = find_max_fft(data_freq);
    freqs_f = sum(data_cnt .* exp(-1i * 2 * pi * max_f / Nfft * [0:Nfft-1]));
%     tmp = [tmp freqs_f];
    if abs(freqs_f) >= thd
        pass_thd_count = pass_thd_count + 1;
    else
        pass_thd_count = 0;
    end
    if pass_thd_count >= 2
        cfo_dec = max_f / Nfft;
        sig_pos = (idx-1) * Nfft + 1;
        break;
    end
end

if sig_pos > 0
    if cfo_dec > 0.5
    cfo_dec = cfo_dec - 1;
    end
    CFO_est = cfo_dec * fsa;
    % CFO compensation
    sig_end_idx = min(sig_pos + 0.52*fsa - 1, length(data_in));
    nn = 0 : (sig_end_idx - sig_pos);
    data_out = data_in(sig_pos : sig_end_idx) .*  exp(-1i * 2 * pi * cfo_dec .* nn);
else
    CFO_est = 0;
    data_out = data_in;
end

% if length(data_coarsed) < n
%     data_coarsed = [data_in(1:sig_start_idx-1) data_coarsed data_in(sig_end_idx+1:n)];
% end
% figure
% plot(abs(tmp))
end

