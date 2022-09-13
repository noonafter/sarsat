clearvars, clc, close all
%% Parameters
fsa = 76800;                     % Sampling clock frequency
fsy = 9600;                      % Symbol frequency 
Ksa = fsa / fsy;                 % Oversampling factor
chan_delay = 0.01;               % Channel delay in sec
chan_atte  = 0.1;                % Channel attenuation
CFO_init   = 13000;                 % Carrier frequency offset(-38.4KHz~38.4KHz) at signal start???
CFO_scan   = -300;               % Carrier frequency offset scanning(-300Hz~300Hz)
CPO_init   = 0.23*pi;             % Carrier phase offset at signal start
SNR_set    = 6:15;               % Signal-to-Noise Ratio of baseband
loop_num   = 1000;                 % Transmission numbers(One frame per transmission)
FIG_DISP   = 'ON';               % MER figure display switch, ON or OFF
%% System Config
% Frame format
pre_len = 136; % for temp
info_len = 192;
check_len = 56;
% BCH generator polynomial
gx = [1 1 1, 1 1 0, 0 1 1, 0 1 1, 0 0 0, 0 1 1, 0 0 1, 0 1 0, 1 1 1, 0 0 0, 1 0 0, 0 1 0, 0 0 0, 1 1 1, 0 1 0, 0 1 0, 0 1 1, 1 0 0, 0 0 1]; % 7633 0312 7042 0722 341
% Preamble
[pre_i,pre_q] = mseq_gene(38400 + pre_len/2);
pre_cplx = pre_i(end-pre_len/2+1:end) + 1i * pre_q(end-pre_len/2+1:end);
% Tx Filter (RRC)
span = 8;
rfc = 0.4;
lpf_coe = rcosdesign(rfc, span, Ksa);
num_snr = length(SNR_set);
bit_errors = zeros(1,num_snr);
ferrors = zeros(1,num_snr);
fsyncs = zeros(1,num_snr);
% err_cs = [];
% err_fs = [];
tic
for idx_snr = 1:num_snr
SNR = SNR_set(idx_snr);
for idx_loop = 1:loop_num
%% --------------------------- Tx ---------------------------------
% Information input
info_bits = randi([0,1], 1, info_len);
check_bits = poly_mod([info_bits zeros(1,check_len)], gx);
pld = [info_bits check_bits];
% Grouping
pld_i = pld(1:2:end);
pld_q = pld(2:2:end); 
% Framing
tagc = 32;
frame_i = [pld_i(1:tagc) pre_i(end-pre_len/2+1:end) pld_i];
frmae_q = [pld_q(1:tagc) pre_q(end-pre_len/2+1:end) pld_q];
% Mapping
frame_bit_i = [2 * frame_i - 1 zeros(1, floor(length(lpf_coe)/2) )];
frame_bit_q = [2 * frmae_q - 1 zeros(1, floor(length(lpf_coe)/2) )];
% Pulse Shaping
sample_i = sqrt(2) / 2 * filter(lpf_coe, 1, sqrt(Ksa) * upsample(frame_bit_i, Ksa));
sample_q = sqrt(2) / 2 * filter(lpf_coe, 1, sqrt(Ksa) * upsample(frame_bit_q, Ksa));
tx_cplx = sample_i + 1i * sample_q;
pre_cplx_up  = tx_cplx(ceil(length(lpf_coe)/2)+tagc*Ksa + [0:pre_len/2*Ksa - 1]);
fsa_pre = fsa; 
%% -------------------------- Channel -----------------------------
% Add CFO
tt = 1 / fsa * [0:length(sample_i)-1];
f_offset = CFO_init + CFO_scan * tt;
tx_cfo = tx_cplx * exp(1i * CPO_init) .*  exp(1i * 2 * pi * f_offset .* tt);
% Channel Delay
tx_delay = [zeros(1,floor(chan_delay*fsa)) tx_cfo];
% AWGN
snr_allband = SNR - 10 * log10(fsa/2 / 24000);
noise_power = 0 - snr_allband;
noise = wgn(1, length(tx_delay), noise_power, 'complex');
% rx_cplx = chan_atte * awgn(tx_delay, snr_allband, 'measured');
rx_cplx = chan_atte * (noise + tx_delay);
%% --------------------------- Rx ----------------------------------
% AGC
agc_coe = 1/6;
desired_mag = 1;
gmin = 2;
gmax = 30;
[data_agced,gain_set] = agc_loop(rx_cplx, agc_coe, desired_mag, gmin, gmax);
% Coarse sync & CFO Estimate 
delay = Ksa;
wlen = length(pre_cplx_up)-delay;
thd = 10^(-1/10) / (1 + 10^(-1/10));
Nfft = 512;
[data_coarsed, sig_start_idx, CFO_est] = cfo_coarse_t018(data_agced, fsa, ...
                                                     pre_cplx_up, fsa_pre, ...
                                                     delay, wlen, thd, Nfft);
if sig_start_idx < 0
    continue;
end
% err_cs = [err_cs CFO_est - CFO_init];
% err_fs = [err_fs sig_start_idx - chan_delay * fsa - ceil(length(lpf_coe)/2) - tagc * Ksa];
% Match Filter
data_matched = filter(lpf_coe, 1, data_coarsed);
% Frame sync & Symbol sync 
Ksa_pre = fsa / fsa_pre;
frame_sync_cplx = pre_cplx_up;
fsync_win_len = length(frame_sync_cplx);
search_len = length(lpf_coe);
fcor = zeros(1,search_len);
for idx_fsync = 1:search_len
    data_fsync_cnt = data_matched(idx_fsync + Ksa_pre * [0:fsync_win_len - 1]);
    fcor(idx_fsync) = data_fsync_cnt * frame_sync_cplx';
end
[~,fmax_idx_relative] = max(abs(fcor));
fmax_idx = fmax_idx_relative ; % where frame_sync_seq begin as to data_matched
pld_idx = fmax_idx + length(frame_sync_cplx) * Ksa_pre;
if pld_idx + (info_len + check_len)/2 * Ksa > length(data_matched)
    continue;
end
% CPO est
nums_cpo = 32;
data_cpo = data_matched( pld_idx - nums_cpo * Ksa + Ksa * [0:nums_cpo-1]);
data_pre = pre_cplx(pre_len/2 - nums_cpo + 1 : end);
cplx_cpo = data_cpo * data_pre';
cpo_est = angle(cplx_cpo);
% Demodulation
demo_info = oqpsk_demod_bst(data_matched, pld_idx, cpo_est);
% BCH decode
dec_info = bchdec_bst(demo_info);
% Record
exp_pos = chan_delay * fsa + (tagc + pre_len/2) * Ksa + length(lpf_coe);
if chan_delay * fsa + tagc * Ksa + pld_idx + floor(length(lpf_coe)/2) == exp_pos
    fsyncs(idx_snr) = fsyncs(idx_snr) + 1;
end
bit_error = sum(abs(dec_info - info_bits));
if bit_error > 0
    bit_errors(idx_snr) = bit_errors(idx_snr) + bit_error;
    ferrors(idx_snr) = ferrors(idx_snr) + 1;
end
end
end
%% MER & Plot
if strcmp(FIG_DISP,'ON') == true
figure('Name', 'AGC增益曲线');
plot((0:length(gain_set)-1)/fsa,gain_set)
figure('Name','OQPSK接收端输入信号星座图')
scatter(real(rx_cplx(1:end-1)),imag(rx_cplx(1:end-1)))
figure('Name','帧同步相关峰搜索曲线')
plot(abs(fcor))
figure('Name','OQPSK匹配滤波前时域信号')
plot((0:length(data_coarsed)-1)/fsa,real(data_coarsed));hold on
plot((0:length(data_coarsed)-1)/fsa,imag(data_coarsed));
figure('Name','OQPSK匹配滤波后时域信号')
plot((0:length(data_matched)-1)/fsa,real(data_matched));hold on
plot((0:length(data_matched)-1)/fsa,imag(data_matched));
figure('Name','OQPSK匹配滤波前信号功率谱密度')
pwelch(data_coarsed);
figure('Name','OQPSK匹配滤波后信号功率谱密度')
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
toc