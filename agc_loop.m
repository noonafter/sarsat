function [data_out,gain_set] = agc_loop(data_in,agc_coe,desired_mag,min_gain,max_gain)
%AGC_LOOP 此处显示有关此函数的摘要
%   此处显示详细说明

n = length(data_in);
data_out = zeros(1,n);
gain_set = zeros(1,n);
gain_cnt = min_gain;
for idx = 1:n
    data_out(idx) = data_in(idx) * gain_cnt;
    gain_cnt = gain_cnt + agc_coe * (desired_mag - abs(data_out(idx)));

    if gain_cnt > max_gain
        gain_cnt = max_gain;
    end
    
    if gain_cnt < min_gain
        gain_cnt = min_gain;
    end
    
    gain_set(idx) = gain_cnt;
end

end

