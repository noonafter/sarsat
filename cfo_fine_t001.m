function [data_out, CFO_fine_est, CFO_scan_est] = cfo_fine_t001(data_in, fsa, Ksa, ucd_search)
%CFO_FINE2_T001 此处显示有关此函数的摘要
%   此处显示详细说明
sig_pos = 1; % data_in only contain valid data

% CFO estimation
n = length(data_in);
search_len = ucd_search * fsa;
win_width = 4 * Ksa;
win_interval = Ksa;
slide_step = Ksa/2;

data_diffed = delay_diff(data_in, win_interval);
data_diffed = data_diffed(win_interval+1:end);
max_pos = search_len - win_width - win_interval + 1;
stop_idx = floor((max_pos - sig_pos)/slide_step) + 1;
cors = zeros(1,stop_idx);
for idx_cor = 1:stop_idx
    data_cnt = data_diffed(sig_pos + (idx_cor-1)*slide_step +  [0:win_width-1]);
    cors(idx_cor) = sum(data_cnt);
end
cfo_set = angle(cors) * fsa / 2 / pi / win_interval;
% phase_cors = angle(cors);
% phase_tmp = phase_cors(1);
% for idx_proc = 2:length(phase_cors)
%     if phase_cors(idx_proc) - phase_tmp > pi
%        phase_cors(idx_proc) = phase_cors(idx_proc) - 2*pi;
%     elseif phase_cors(idx_proc) - phase_tmp < -pi
%         phase_cors(idx_proc) = phase_cors(idx_proc) + 2*pi;
%     end
%     phase_tmp = phase_cors(idx_proc);
% end
% cfo_set = phase_cors * fsa / 2 / pi / win_interval;

rm = length(cfo_set);
% b = regress(cfo_set',[ones(rm,1) (0:rm-1)']);
yk = sum(cfo_set);
kyk = (0:rm-1) * cfo_set';
sca = 2/(rm * (rm-1) * (rm+1));
b = sca * [(rm-1)*(2*rm-1) -3*(rm-1); -3*(rm-1) 6] * [yk;kyk];

k = b(2) * fsa / slide_step / 2;
f0 = b(1) - k * Ksa / fsa;
CFO_fine_est = f0;
CFO_scan_est = k;

% CFO compensation
com_start_idx = sig_pos + search_len;
com_end_idx = min(com_start_idx + (0.52-0.1)*fsa - 1, n);
tt = ucd_search + 1/fsa * [0 : com_end_idx - com_start_idx];
cfo_fine = f0 + k * tt;
data_out = data_in(com_start_idx : com_end_idx) .*  exp(-1i * 2 * pi * cfo_fine .* tt);

% if length(data_out) < n
%     data_out = [data_in(1:com_start_idx-1) data_out data_in(com_end_idx+1:n)];
% end
end

