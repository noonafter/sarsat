clearvars, clc, close all
%% Parameters
fsa        = 153600;             % Sampling clock frequency
fsy        = 38400;              % Symbol frequency (Chip rate)
Ksa        = fsa / fsy;          % Oversampling factor
Kss        = 256;                % Spread spectrum factor of I/Q path
chan_delay = 0.01;               % Channel delay in sec
chan_atte  = 0.1;                % Channel attenuation
CFO_init   = 13000;              % Carrier frequency offset(-38.4kHz~38.4KHz) at signal start
CFO_scan   = -300;               % Carrier frequency offset scanning(-300Hz~300Hz)
CPO_init   = 0.3*pi;             % Carrier phase offset at signal start
SNR_set    = -6:3;               % Signal-to-Noise Ratio of baseband
loop_num   = 1000;                 % Transmission numbers(One frame per transmission)
FIG_DISP   = 'ON';               % MER figure display switch, ON or OFF
%% System Config
% Frame format
pre_len = 50;
info_len = 202;
check_len = 48;
pre_chip_len = pre_len * Kss / 2;
total_len = pre_len + info_len + check_len;
% DSSS
[mseq_i,mseq_q] = mseq_gene(fsy * 1);
sseq_i = 2 * mseq_i - 1;
sseq_q = 2 * mseq_q - 1;
sseq_i_tmp = sum([sseq_i; sseq_i(2:end) -sseq_i(end)])/2;
sseq_q_tmp = sum([sseq_q; sseq_q(2:end) -sseq_q(end)])/2;
sseq_upi = [reshape([sseq_i; sseq_i_tmp], 1, []) 0];
sseq_upq = [0 reshape([sseq_q; sseq_q_tmp], 1, [])];
sseq_up_cplx = sseq_upi + 1i * sseq_upq;
pre_ups_cplx = sseq_up_cplx(1 : 2 * pre_chip_len);
pld_ups_cplx = sseq_up_cplx(2 * pre_chip_len + 1 : end);
fsa_pre = fsy * 2; 
% BCH generator polynomial
gx = [1 1 1 0, 0 0 1, 1 1 1, 1 1 0, 1 0 1, 1 1 0, 0 0 0, 1 0 1, 1 1 0, 1 1 1, 1 1 0, 0 1 1, 1 1 0, 0 1 0, 0 1 0, 1 1 1];
% Tx Filter (RRC)
span = 8;
lpf_coe = rcosdesign(0.4,span,Ksa);
num_snr = length(SNR_set);
bit_errors = zeros(1,num_snr);
ferrors = zeros(1,num_snr);
fsyncs = zeros(1,num_snr);
err_cs = [];
err_fs = [];
tic
for idx_snr = 1:num_snr
SNR = SNR_set(idx_snr);
for idx_loop = 1:loop_num
%% --------------------------- Tx ---------------------------------
% Information input
info_bits = randi([0,1],1,info_len);
% -------------------------------FOR TEST ONLY-------------------------- %
% info_bits = [0 0 0 0 0 0 0 0 1 1 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 1 0 1 0 0 ...
% 1 1 0 0 1 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 1 0 1 1 0 0 0 0 1 ...
% 1 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 0 0 0 0 0 0 ...
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
% 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 ...
% 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 ...
% 1 0 0 1 0 1 1 0 0 0];
% -------------------------------FOR TEST ONLY-------------------------- %
% BCH coding
check_bits = poly_mod([info_bits zeros(1,check_len)], gx);
% Framing
frame_all = [zeros(1,pre_len) info_bits check_bits];
% Grouping
frame_i = frame_all(1:2:end);
frame_q = frame_all(2:2:end);
% Spread spectrum
frame_i_tmp = reshape(repmat(frame_i,Kss,1),1,[]);
frame_q_tmp = reshape(repmat(frame_q,Kss,1),1,[]);
frame_dsss_i = mod(mseq_i + frame_i_tmp,2);
frame_dsss_q = mod(mseq_q + frame_q_tmp,2);
% Mapping
frame_bit_i = [2 * frame_dsss_i - 1 zeros(1, length(lpf_coe) )];
frame_bit_q = [2 * frame_dsss_q - 1 zeros(1, length(lpf_coe) )];
% Pulse shaping
sample_i_tmp = sqrt(2) / 2 * filter(lpf_coe, 1, sqrt(Ksa) * upsample(frame_bit_i,Ksa));
sample_q_tmp = sqrt(2) / 2 * filter(lpf_coe, 1, sqrt(Ksa) * upsample(frame_bit_q,Ksa));
% OQPSK
sample_i = [sample_i_tmp   zeros(1,Ksa/2)];
sample_q = [zeros(1,Ksa/2)   sample_q_tmp]; % check energy and spectrum
tx_cplx = sample_i + 1i * sample_q;
%% -------------------------- Channel -----------------------------
% Add CFO
tt = 1 / fsa * [0:length(sample_i)-1];
f_offset = CFO_init + CFO_scan * tt;
tx_cfo = tx_cplx * exp(1i * CPO_init) .*  exp(1i * 2 * pi * f_offset .* tt);
% Channel Delay
tx_delay = [zeros(1,floor(chan_delay*fsa)) tx_cfo];
% AWGN
snr_allband = SNR - 10 * log10(fsa/2 / 40000);
rx_cplx = chan_atte * awgn(tx_delay, snr_allband, 'measured');
%% --------------------------- Rx ----------------------------------
% AGC
agc_coe = 1/32;
desired_mag = 1;
gmin = 1;
gmax = 30;
[data_agced,gain_set] = agc_loop(rx_cplx, agc_coe, desired_mag, gmin, gmax);
% Coarse sync & CFO Estimate
wst = 512;
pre_coarse = pre_ups_cplx(wst+1:end);
delay = 4;
wlen = 2048;
thd = 10^(-9/10) / (1 + 10^(-9/10));
Nfft = 512;
[data_coarsed, sig_start_idx, CFO_coarse_est] = cfo_coarse_t018(data_agced, fsa, ...
                                                            pre_coarse, fsa_pre, ...
                                                            delay, wlen, thd, Nfft);
% err_cs = [err_cs CFO_coarse_est-(CFO_init+CFO_scan*6.67/1000)];
% err_fs = [err_fs sig_start_idx-(chan_delay*fsa+wst*2+ceil(length(lpf_coe)/2))];
if sig_start_idx < 0 || sig_start_idx + (1-0.166)*fsa - 1 > length(data_agced) 
    continue;
end
% Match Filter
data_matched = filter(lpf_coe, 1, data_coarsed);
data_matched = data_matched(ceil(length(lpf_coe)/2):end);
% Fine CFO Estimate
wst_tail = 1 * Kss*fsa_pre/fsy;
pre_fine = pre_ups_cplx(wst + 1 : end - wst_tail);
delay = 512;
wlen = 512;
step = 256;
[data_fined, CFO_fine_est, CFO_scan_est] = cfo_fine_t018(data_matched, fsa, ...
                                                    pre_fine, fsa_pre, ...
                                                    delay, wlen, step);
% err_cs = [err_cs CFO_fine_est-(CFO_init+CFO_scan*6.67/1000 - CFO_coarse_est)];
% err_fs = [err_fs CFO_scan_est - CFO_scan];
% Frame sync & Symbol sync 
pre_sync = pre_ups_cplx(end - wst_tail + 1:end);
fsync_win_len = length(pre_sync);
search_len = 4 * length(lpf_coe);
fcor = zeros(1,search_len);
for idx_fsync = 1 : search_len
    data_fsync_cnt = data_fined(idx_fsync + Ksa/2 * [0:fsync_win_len - 1]);
    fcor(idx_fsync) = data_fsync_cnt * pre_sync';
end
[~,fmax_idx] = max(abs(fcor));
pld_idx = fmax_idx + wst_tail * fsa_pre/fsy;
if pld_idx + (info_len + check_len)/2 * Kss * Ksa > length(data_fined)
    continue;
end
% CPO est
cplx_cpo = fcor(fmax_idx);
cpo_est = cplx_cpo / abs(cplx_cpo);
% Demodulation
demo_info = oqpsk_demod_t018(data_fined, pld_idx, cpo_est, pld_ups_cplx);
% BCH decode
dec_info = bchdec_t018(demo_info);
% Record
exp_pos = chan_delay * fsa + pre_len * Kss/2 * Ksa + length(lpf_coe); % pld expectation position
real_pos = sig_start_idx-1 + floor(length(lpf_coe)/2) + (length(data_matched)-length(data_fined)) + pld_idx;
if  real_pos >= exp_pos - 2 && real_pos <= exp_pos + 2
    fsyncs(idx_snr) = fsyncs(idx_snr) + 1;
end
bit_error = sum(abs(dec_info - info_bits));
if bit_error > 0
    bit_errors(idx_snr) = bit_errors(idx_snr) + bit_error;
    ferrors(idx_snr) = ferrors(idx_snr) + 1;
end
end
end
toc
%% MER & Plot
if strcmp(FIG_DISP,'ON') == true
figure('Name', 'AGC增益曲线');
plot((0:length(gain_set)-1)/fsa,gain_set)
figure('Name','T018接收端输入信号星座图')
scatter(real(rx_cplx(1:end-1)),imag(rx_cplx(1:end-1)))
figure('Name','帧同步相关峰搜索曲线')
plot(abs(fcor))
figure('Name','T018匹配滤波前时域信号')
plot((0:length(data_coarsed)-1)/fsa,real(data_coarsed));hold on
plot((0:length(data_coarsed)-1)/fsa,imag(data_coarsed));
figure('Name','T018匹配滤波后时域信号')
plot((0:length(data_matched)-1)/fsa,real(data_matched));hold on
plot((0:length(data_matched)-1)/fsa,imag(data_matched));
figure('Name','T018匹配滤波前信号功率谱密度')
pwelch(data_coarsed);
figure('Name','T018匹配滤波后信号功率谱密度')
pwelch(data_matched);
figure('Name', '帧同步成功率');
plot(SNR_set, fsyncs / loop_num)
xlabel('SNR(dB)');
ylabel('SSR');
figure('Name', '误帧率');
plot(SNR_set, ferrors / loop_num)
% plot(SNR_set, (ferrors + (loop_num - fsyncs))/ loop_num)
xlabel('SNR(dB)');
ylabel('FER');
figure('Name', '误比特率');
plot(SNR_set, bit_errors / (loop_num * length(info_bits))) 
% plot(SNR_set, (bit_errors + (loop_num - fsyncs)* length(info_bits)) / (loop_num * length(info_bits)))
xlabel('SNR(dB)');
ylabel('BER');
end
