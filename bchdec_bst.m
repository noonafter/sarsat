function [decode] = bchdec_bst(code)
%BCHDEC_T018 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
rx_code = gf([zeros(1,7) code]);
deco_rx = bchdec(rx_code, 255, 199);
decode = double(deco_rx.x(8:end));
end

