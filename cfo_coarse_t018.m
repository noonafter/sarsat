function [data_out, sig_pos, CFO_est] = cfo_coarse_t018(data_in, fsa, preamble, fsa_pre, delay, wlen, thd, Nfft)
%CFO_COARSE3_T018 此处显示有关此函数的摘要
%   此处显示详细说明
Ksa_pre = fsa / fsa_pre;
sig_pos = -1;
cor_tmp = 0;
% tmp = [];
% Difference
data_diffed = delay_diff(data_in, delay * Ksa_pre);
pmb_diffed  = delay_diff(preamble, delay);
pmb_valid = pmb_diffed(delay+(1:wlen));
nrom_fac = sqrt(pmb_valid * pmb_valid');
pmb_valid = pmb_valid / nrom_fac;
idx_stop = length(data_diffed) - Ksa_pre * wlen + 1;
% Coarse Sync
for idx = 1:idx_stop
    data_cnt = data_diffed(idx + Ksa_pre * [0:wlen - 1]);
    cor = data_cnt * pmb_valid'/ sqrt(data_cnt*data_cnt');
%     tmp = [tmp cor];
    if abs(cor) > thd
        if abs(cor) > abs(cor_tmp)
            cor_tmp = cor;
            sig_pos = idx;
            continue;
        else
            break;
        end
    elseif sig_pos > 0
        break;
    end
end

if sig_pos > 0
    % CFO estimation
    sig_pos = sig_pos - delay * Ksa_pre;
%     CFO_est = angle(cor_tmp) / 2 / pi / (delay * Ksa/2) * fsa;
    tone_demo = data_in(sig_pos + Ksa_pre * [0:Nfft-1]) .* conj(preamble(1 + [0:Nfft-1]));
    data_freq = fft(tone_demo, Nfft);
    max_f = find_max_fft(data_freq);
    cfo_severe_est = max_f / Nfft;
    if cfo_severe_est > 0.5
        cfo_severe_est = cfo_severe_est - 1;
    end
    CFO_est = cfo_severe_est * fsa/Ksa_pre ;
    % CFO compensation
    sig_end_idx = min(sig_pos + 1*fsa - 1, length(data_in));
    nn = 0 : (sig_end_idx - sig_pos);
    data_out = data_in(sig_pos : sig_end_idx) .*  exp(-1i * 2 * pi * CFO_est .* nn / fsa);
else
    CFO_est = 0;
    data_out = data_in;
end

% figure
% plot(abs(tmp))
end