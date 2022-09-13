function [decode] = bchdec_t018(code)
%BCHDEC_T018 此处显示有关此函数的摘要
%   此处显示详细说明
rx_code = gf([zeros(1,5) code]);
deco_rx = bchdec(rx_code, 255, 207);
decode = double(deco_rx.x(6:end));
end

