function [data_out, CFO_est, CFS_est] = cfo_fine_t018(data_in, fsa, preamble, fsa_pre, delay, wlen, step)
%CFO_FINE2_T018 此处显示有关此函数的摘要
%   此处显示详细说明
Ksa_pre = fsa / fsa_pre; % 2 samples : 1 preamble sym
data_tone = data_in(1 + Ksa_pre * (0:length(preamble)-1)) .* conj(preamble);
data_diffed = delay_diff(data_tone, delay);
data_diffed = data_diffed(delay+1:end);
cors = [];
left_idx = 1;
stop_idx = length(data_diffed) - wlen + 1 ; % the last 256 chip is not used
while left_idx <= stop_idx
    data_cnt = data_diffed(left_idx + [0:wlen-1]);
    cors = [cors sum(data_cnt)];
    left_idx = left_idx + step;
end
phase_cors = angle(cors);
phase_tmp = phase_cors(1);
for idx_proc = 2:length(phase_cors)
    if phase_cors(idx_proc) - phase_tmp > pi
       phase_cors(idx_proc) = phase_cors(idx_proc) - 2*pi;
    elseif phase_cors(idx_proc) - phase_tmp < -pi
        phase_cors(idx_proc) = phase_cors(idx_proc) + 2*pi;
    end
    phase_tmp = phase_cors(idx_proc);
end
cfo_set = phase_cors * fsa_pre / 2 / pi / delay;

% Linear Regression
rm = length(cfo_set);
% b = regress(cfo_set',[ones(rm,1) (0:rm-1)']);
yk = sum(cfo_set);
kyk = (0:rm-1) * cfo_set';
sca = 2/(rm * (rm-1) * (rm+1));
b = sca * [(rm-1)*(2*rm-1) -3*(rm-1); -3*(rm-1) 6] * [yk;kyk];

k = b(2) * fsa_pre / step / 2;
f0 = b(1)  - k * delay / fsa_pre;
CFO_est = f0;
CFS_est = k;
% CFO compensation
com_start_idx = (length(preamble)-16) * Ksa_pre; % compensate another 16, to see complete cor peak
com_end_idx = min(com_start_idx + (1-0.1)*fsa - 1, length(data_in));
tt = 1/fsa * [com_start_idx : com_end_idx];
cfo_fine = f0 + k * tt;
data_out = data_in(com_start_idx : com_end_idx) .*  exp(-1i * 2 * pi * cfo_fine .* tt);
end

