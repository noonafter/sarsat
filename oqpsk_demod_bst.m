function [info_out] = oqpsk_demod_bst(data_in, pld_idx, cpo_init)
%OQPSK_DEMOD_BST 此处显示有关此函数的摘要
%   此处显示详细说明
% System Config
info_len = 192;
check_len = 56;
Ksa = 8;
% Demodulation
Kp = 1;
cpo_cnt = cpo_init;
demod_bits = zeros(2,(info_len + check_len)/2);
for idx_sym = 1:(info_len + check_len)/2
    data_sym_cnt = data_in(pld_idx + (idx_sym-1) * Ksa);
    data_cpoed = data_sym_cnt * exp(-1i * cpo_cnt);
    slice_i = real(data_cpoed);
    slice_q = imag(data_cpoed);
    if slice_i > 0
        demod_bits(1,idx_sym) = 1;
        data_tx_i = 1;
    else
        demod_bits(1,idx_sym) = 0;
        data_tx_i = -1;
    end
    if slice_q > 0
        demod_bits(2,idx_sym) = 1;
        data_tx_q = 1;
    else
        demod_bits(2,idx_sym) = 0;
        data_tx_q = -1;
    end
    ped = angle(data_cpoed * (data_tx_i - 1i * data_tx_q));
    cpo_cnt = cpo_cnt + Kp * ped;
end
info_out = reshape(demod_bits,1,[]);
end

