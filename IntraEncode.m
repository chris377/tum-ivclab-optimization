function dst = IntraEncode(image, qScale)
%  Function Name : IntraEncode.m
%  Input         : image (Original RGB Image)
%                  qScale(quantization scale)
%  Output        : dst   (sequences after zero-run encoding, 1xN)
 
    zz_all = zeros(1,numel(image));
    idx_zz_all = 1;
    % DCT
    coeff = blockproc(image,[8 8],@(block_struct) DCT8x8(block_struct.data));
    % Quantization
    quant = blockproc(coeff,[8 8],@(block_struct) Quant8x8(block_struct.data, qScale));
    
    % zigzag scan
    for j = 1: 8: size(image,2)
        for i = 1: 8: size(image,1)
            loc_imageBlock = quant(i:i+7,j:j+7,:);
            zz = ZigZag8x8(loc_imageBlock );              
            zz_all(1,idx_zz_all:idx_zz_all+191) = zz(:)'; %  transforn into vector 1 X 64*3 = 1 X 192
            idx_zz_all = idx_zz_all + 192;
        end
    end
    dst = ZeroRunEnc_EoB(zz_all(:)', 1000);
end
