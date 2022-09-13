function [dec_info] = bchdec_t001(soft_info, dec_type)
%BCHDEC_T001 此处显示有关此函数的摘要
%   此处显示详细说明
dec_type = 0; % only support hard decode for this version

if length(soft_info) < 144 % prevent idx overflow
    soft_info = [soft_info zeros(1,144 - length(soft_info))];
end

if dec_type  == 0
    demo_info = soft_info > 0; % hard slicer
    if demo_info(25) == 0 % short msg
        rxcode_pdf1 = gf([zeros(1,45) demo_info(25:106)]);
        rxcode_npdf = demo_info(107:112);
        decoded_pdf1 = bchdec(rxcode_pdf1, 127, 106);
        dec_info = [decoded_pdf1.x(46:106) rxcode_npdf];
    else
        rxcode_pdf1 = gf([zeros(1,45) demo_info(25:106)]);
        rxcode_pdf2 = gf([zeros(1,25) demo_info(107:144)]);
        decoded_pdf1 = bchdec(rxcode_pdf1, 127, 106);
        decoded_pdf2 = bchdec(rxcode_pdf2, 63, 51);
        dec_info = [decoded_pdf1.x(46:106) decoded_pdf2.x(26:51)];
    end
else
    % wait for bch soft decode...
    
end
dec_info = double(dec_info);



end

