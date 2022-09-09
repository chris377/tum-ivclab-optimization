function dst = IntraDecode(image, img_size , qScale)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image)

    % RL Dec
    recZZ = ZeroRunDec_EoB(image,1000);
    loc_dst = zeros(img_size);
    idx = 1;

    % iZZ
    for j = 1:8:img_size(1,2)
        for i = 1:8:img_size(1,1)
            recZZBlock = recZZ(idx:idx+191); % 8*8*3 Block
            rep_recZZBlock =reshape(recZZBlock,[64,3]); % 64*3 Block         
            quant_block = DeZigZag8x8(rep_recZZBlock);   
            loc_dst(i:i+7,j:j+7,:) = quant_block;
            idx = idx + 192;
        end
    end
    
    % deQuan
    dct_block = blockproc(loc_dst,[8 8],@(block_struct) DeQuant8x8(block_struct.data, qScale));
        
    % Idct
    dst  =  blockproc(dct_block,[8 8],@(block_struct) IDCT8x8(block_struct.data));

end