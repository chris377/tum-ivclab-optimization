function block = IDCT8x8(coeff)
%  Function Name : IDCT8x8.m
%  Input         : coeff (DCT Coefficients) 8*8*3
%  Output        : block (original image block) 8*8*3
    S = idct(coeff,[],2);
    T = idct(S,[],1);
    
    % Update value
    block = T;

end