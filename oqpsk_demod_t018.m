function [info_out] = oqpsk_demod_t018(data_in, pld_idx, cpo_init, pld_ups_cplx)
%OQPSK_DEMOD_T018 此处显示有关此函数的摘要
%   此处显示详细说明
% System Config
info_len = 202;
check_len = 48;
Kss = 256;
Ksa = 4;
% Demodulation
cpo_cnt = cpo_init;
demod_bits = zeros(2,(info_len + check_len)/2);
last_q = 0;
% Kp = 1;
% nn = 0:Kss*2;
for idx_sym = 1:(info_len + check_len)/2
    data_sym_cnt = data_in(pld_idx + (idx_sym-1) * Kss * Ksa + Ksa/2 * (0:Kss*2));
    data_cpoed = data_sym_cnt * conj(cpo_cnt);
%     data_cpoed = data_sym_cnt * exp(-1i * cpo_cnt) .* exp(-1i * cfo_cnt * nn / fsa*2);
    data_qpsk_i = real(data_cpoed(1:Kss*2));
    data_qpsk_q = imag(data_cpoed(2:Kss*2+1));
    mseq_cnt_i = real(pld_ups_cplx((idx_sym-1) * Kss*2 + (1:Kss*2)));
    mseq_cnt_q = imag(pld_ups_cplx((idx_sym-1) * Kss*2 + (2:Kss*2+1)));
    slicer_i = data_qpsk_i * mseq_cnt_i';
    slicer_q = data_qpsk_q * mseq_cnt_q';
    if slicer_i > 0
        demod_bits(1,idx_sym) = 0;
        data_tx_i = mseq_cnt_i;
    else
        demod_bits(1,idx_sym) = 1;
        data_tx_i = -mseq_cnt_i;
    end
    if slicer_q > 0
        demod_bits(2,idx_sym) = 0;
        data_tx_q = [last_q mseq_cnt_q(1:end-1)];
        last_q = mseq_cnt_q(end);
    else
        demod_bits(2,idx_sym) = 1;
        data_tx_q = [last_q -mseq_cnt_q(1:end-1)];
        last_q = -mseq_cnt_q(end);
    end
    % update cpo
    data_tx = data_tx_i + 1i * data_tx_q;
    ped = data_cpoed(1:end-1) * data_tx';
    ped_norm = ped / abs(ped);
    cpo_cnt = cpo_cnt * ped_norm;
end
info_out = reshape(demod_bits,1,[]);
end

