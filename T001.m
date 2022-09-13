clearvars, clc, close all
%% Parameters
fsy        = 400;               % Symbol frequency(Bit rate for binary modulation)
fsa        = 38400;             % Sampling clock frequency
Ksa        = fsa / fsy;         % Oversampling factor
CFO_init   = -13000;             % Carrier frequency offset(-19.2kHz~19.2KHz) at signal start
CFO_scan   = +300;              % Carrier frequency offset scanning(-600Hz~600Hz)
CPO_init   = 0.3*pi;            % Carrier phase offset at signal start
msg_type   = 1;                 % 0: short, 1: long
dec_type   = 0;                 % 0: hard decode, 1: soft decode
chan_delay = 0.05;              % Channel delay in sec
chan_atte  = 0.1;               % Channel attenuation(Including Tx, Channel, Rx RF gain)
SNR_set    = 4:13;              % Ratio of signal to noise power
loop_num   = 10;               % Transmission numbers(One frame per transmission)
FIG_DISP   = 'ON';              % MER figure display switch, ON or OFF
%% System Config
ucd = 0.16;
gx1 = [1 0 0 1 1 0 1 1 0 1 1 0 0 1 1 1 1 0 0 0 1 1];
gx2 = [1 0 1 0 1 0 0 1 1 1 0 0 1];
bit_seq = ones(1,15);
frame_seq = [0 0 0 1 0 1 1 1 1];
mag_h = tan(0.99);mag_t = tan(1.1);
pulse_coe_i = [(0:4)/5 ones(1,44) (4:-1:0)/5];
pulse_norm_i =  pulse_coe_i ./ sqrt(sum(abs(pulse_coe_i).^2));
pulse_coe_q = [(0:4)*mag_h/5 interp1([6 27 28 49],[mag_h mag_t mag_t mag_h],6:49,'PCHIP') (4:-1:0)*mag_h/5];
pulse_norm_q =  pulse_coe_q ./ sqrt(sum(abs(pulse_coe_q).^2));
matched_coe_i = flip(pulse_norm_i);
matched_coe_q = flip(pulse_norm_q);
frame_sync_seq = [bit_seq frame_seq];
frame_sync_tmp = 2 * frame_sync_seq - 1;
frame_sync_tmp2 = reshape([frame_sync_tmp; -frame_sync_tmp],1,[]);
frame_sync_cplx = exp(1i * 1.1 .* frame_sync_tmp2);
% MER
snr_num = length(SNR_set);
fsync_set = zeros(1,snr_num);
ferr_set = zeros(1,snr_num);
berr_set = zeros(1,snr_num);
% err_cs = [];err_cf = [];err_ccs = [];
tic
for idx_snr = 1:snr_num
SNR = SNR_set(idx_snr);
for idx_loop = 1:loop_num
%% --------------------------- Tx ------------------------------------
% Information Bits
pdf1 = [msg_type randi([0,1],1,60)];
% -------------------FOR TEST ONLY------------------- %
% pdf1 = [0 1 0 1 0 1 1 0 1 1 1 0 0 1 1 0 1 zeros(1,8) 1 zeros(1,12) 1 0 0 0 1 zeros(1,7) 1 zeros(1,9) 1];
% -------------------FOR TEST ONLY------------------- %
if msg_type == 0  
    npdf = randi([0,1],1,6);
% BCH Coding
    check_pdf1 = poly_mod([pdf1 zeros(1,length(gx1)-1)], gx1);
% Framing
    info_bits = [pdf1 npdf];
    msg = [bit_seq, frame_seq, pdf1, check_pdf1, npdf];
else  
    pdf2 = randi([0,1],1,26);
% -------------------FOR TEST ONLY------------------- %
%     pdf2 = [1 0 0 1 0 1 0 1 1 1 zeros(1,11) 1  0 1 1 1];
% -------------------FOR TEST ONLY------------------- %
% BCH Coding
    check_pdf1 = poly_mod([pdf1 zeros(1,length(gx1)-1)], gx1);
    check_pdf2 = poly_mod([pdf2 zeros(1,length(gx2)-1)], gx2);
% Framing
    info_bits = [pdf1 pdf2];
    msg = [bit_seq, frame_seq, pdf1, check_pdf1, pdf2, check_pdf2];
end
% Biphase Code L & Unmodulated Carrier
frame_bit_tmp = 2 * msg - 1;
frame_bit_q = reshape([frame_bit_tmp; -frame_bit_tmp],1,[]);
frame_bit_i = [ ones(1,ucd*fsy*2)/cos(1.1)  ones(size(frame_bit_q)) zeros(1,1)];
frame_bit_q = [zeros(1,ucd*fsy*2)           frame_bit_q             zeros(1,1)];
% Pulse Shaping
sample_i = cos(1.1) * filter(pulse_norm_i, 1, sqrt(Ksa/2) * upsample(frame_bit_i, Ksa/2));
sample_q = sin(1.1) * filter(pulse_norm_q, 1, sqrt(Ksa/2) * upsample(frame_bit_q, Ksa/2));
%% --------------------------- Channel -------------------------------
% Add CFO
t_bias = 0;
tt = 1 / fsa * [t_bias + 0:t_bias + length(sample_i)-1];
f_offset = CFO_init + CFO_scan * tt;
tx_cplx = ( sample_i + 1i * sample_q) * exp(1i * CPO_init) .*  exp(1i * 2 * pi * f_offset .* tt);
% Channel Delay
tx_delay = [zeros(1,floor(chan_delay*fsa)) tx_cplx];
% AWGN
snr_all_band = SNR - 10 * log10(fsa/2 / 3000);
data_channeled = chan_atte * awgn(tx_delay, snr_all_band, 'measured');
%% --------------------------- Rx ------------------------------------
% % AGC
agc_coe = 1/32;
desired_mag = 1;
gmin = 1;
gmax = 30;
[data_agced,gain_set] = agc_loop(data_channeled,agc_coe,desired_mag,gmin,gmax);
% Coarse sync & CFO Estimate
Nfft = 512;
thd = Nfft / 2 * 0.5;
% thd = 0.25;
[data_coarsed, sig_start_idx, CFO_coarse_est] = cfo_coarse_t001(data_agced, fsa, Nfft, thd);
% err_ccs = [err_ccs  CFO_coarse_est - (CFO_init + 1 * Nfft/fsa*CFO_scan)];
if sig_start_idx < 0 || sig_start_idx + (0.28+0.08*msg_type) * fsa - 1 > length(data_agced)
    continue;
end
% Fine CFO Estimate
ucd_search = ucd - 2 * Nfft/fsa;
[data_fined, CFO_fine_est, CFO_scan_est] = cfo_fine_t001(data_coarsed, fsa, Ksa, ucd_search);
% err_cs = [err_cs  CFO_fine_est - (CFO_init + 1 * Nfft/fsa*CFO_scan - CFO_coarse_est)];
% err_cf = [err_cf CFO_scan_est - CFO_scan];
% Match Filter
data_matched = filter(matched_coe_q, 1, data_fined);
% Frame Sync & Symbol Sync
fsync_win_len = length(frame_sync_cplx);
search_len  = 4*Nfft + length(pulse_norm_q);
fcor = zeros(1,search_len);
for idx_fsync = 1:search_len
    data_fsync_cnt = data_matched(idx_fsync + Ksa/2 * [0:fsync_win_len - 1]);
    fcor(idx_fsync) = (data_fsync_cnt) * (frame_sync_cplx)';
end
[~,fmax_idx_relative] = max(abs(fcor));
fmax_idx = fmax_idx_relative; % where frame_sync_seq begin as to data_matched
% Demodulation
data_sym_synced = downsample(data_matched(fmax_idx:end),Ksa/2); % end!!!
data_sym_synced_p = downsample(data_sym_synced,2);
data_sym_synced_n = downsample(data_sym_synced,2,1);
data_sym_synced_n = [data_sym_synced_n zeros(1,length(data_sym_synced_p) - length(data_sym_synced_n))];
data_diff = data_sym_synced_p .* conj(data_sym_synced_n);
demo_soft_info = imag(data_diff);
% BCH Decoder
dec_info = bchdec_t001(demo_soft_info, dec_type); % only suppot hard decode for now
% Record
if sig_start_idx + ucd_search * fsa + fmax_idx > (ucd+chan_delay)*fsa + 2 && ...
   sig_start_idx + ucd_search * fsa + fmax_idx < (ucd+chan_delay)*fsa + 106
    fsync_set(idx_snr)  = fsync_set(idx_snr) + 1;
    if length(dec_info) == length(info_bits)
        berr_set(idx_snr) = berr_set(idx_snr) + sum(abs(dec_info - info_bits));
    else
        berr_set(idx_snr) = berr_set(idx_snr) + length(info_bits);
    end
    if berr_set(idx_snr) > 0
        ferr_set(idx_snr) = ferr_set(idx_snr) + 1;
    end
end

end
end
toc
%% MER & Plot
% fprintf("\nMeasured MER\n");
% fprintf("\nSend: %d frames, Received: %d frames.\n ",  loop_num , fsync_count);
% fprintf("\nSync Success Ratio: %.2f%%\n", (fsync_count / loop_num)*100);
% fprintf("\nBit Error Ratio: %.2f%%\n", (err_count / loop_num/ length(info_bits))*100);
% fprintf("\nCFO Estimation: %.2f\n", CFO_coarse_est + CFO_fine_est);

if strcmp(FIG_DISP,'ON') == true
figure('Name', 'AGC增益曲线');
plot((0:length(gain_set)-1)/fsa,gain_set)
figure('Name','T001接收端输入信号星座图')
scatter(real(data_channeled(1:end-1)),imag(data_channeled(1:end-1)))
figure('Name','T001符号同步后信号星座图')
scatter(real(data_sym_synced(1:end-1)),imag(data_sym_synced(1:end-1)))
axis([-8,8, -8,8])
figure('Name','帧同步相关峰搜索曲线')
plot(abs(fcor))
figure('Name','T001匹配滤波前时域信号')
plot((0:length(data_fined)-1)/fsa,real(data_fined));hold on
plot((0:length(data_fined)-1)/fsa,imag(data_fined));
figure('Name','T001匹配滤波后时域信号')
plot((0:length(data_matched)-1)/fsa,real(data_matched));hold on
plot((0:length(data_matched)-1)/fsa,imag(data_matched));
figure('Name','T001匹配滤波前信号功率谱密度')
pwelch(data_fined);
figure('Name','T001匹配滤波后信号功率谱密度')
pwelch(data_matched);
figure('Name', '帧同步成功率');
plot(SNR_set, fsync_set / loop_num)
xlabel('SNR(dB)');
ylabel('SSR');
figure('Name', '误帧率');
plot(SNR_set, ferr_set / loop_num)
% plot(SNR_set, (ferr_set + (loop_num - fsync_set))/ loop_num)
xlabel('SNR(dB)');
ylabel('FER');
figure('Name', '误比特数');
plot(SNR_set, berr_set / (loop_num * length(info_bits)))
% plot(SNR_set, (berr_set + (loop_num - fsync_set)* length(info_bits)) / (loop_num * length(info_bits)))
xlabel('SNR(dB)');
ylabel('BER');
end

