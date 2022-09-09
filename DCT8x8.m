function coeff = DCT8x8(block)
%  Input         : block    (Original Image block, 8x8x3)
%
%  Output        : coeff    (DCT coefficients after transformation, 8x8x3)
    Q = dct(block,[],1); % along the rows 
    R = dct(Q,[],2);     % along the columns
    
    % Update value
    coeff = R;
    
end