function [recon_image,dst,pred_coeff] = IntraEncode_new(image, qScale)
    [x,y,z] = size(image);
    zz_all = [];
    pred_coeff = [];
    result = zeros(size(image));

    %case1: leftuppercorner
    block1 = image(1:8,1:8,:);
    mode = 1;
    [result_block,zz,cof] = intra_predictor(block1,mode,qScale,0,0,0);
    result(1:8,1:8,:) = result_block;
    zz_all = [zz_all,zz(:)'];
    %pred_coeff = [pred_coeff,cof];


    %case2: first row
    for i = 2:y/8
        block = image(1:8,(i-1)*8+1:i*8,:);
        mode = 2;
        prev_hori = result(1:8,(i-1)*8,:);
        [result_block,zz,cof] = intra_predictor(block,mode,qScale,prev_hori,0,0,0);
        result(1:8,(i-1)*8+1:i*8,:) = result_block;
        zz_all = [zz_all,zz(:)'];
        pred_coeff = [pred_coeff,cof];
    end
    
    %case3: otherwise
    for i = 2:x/8
        for j = 1:y/8
            if j == 1
                block = image((i-1)*8+1:i*8,1:8,:);
                mode = 3;
                prev_verti = result((i-1)*8,1:8,:);
                prev_verti2 = result((i-1)*8,j*8+1:(j+1)*8,:);
                [result_block,zz,cof] = intra_predictor(block,mode,qScale,0,prev_verti,prev_verti2,0);
                result((i-1)*8+1:i*8,1:8,:) = result_block;
                zz_all = [zz_all,zz(:)'];
                pred_coeff = [pred_coeff,cof];
            elseif j == y/8
                block = image((i-1)*8+1:i*8,(j-1)*8+1:j*8,:);
                prev_verti = result((i-1)*8,(j-1)*8+1:j*8,:);
                prev_hori = result((i-1)*8+1:i*8,(j-1)*8,:);
                diagonalvalue = result((i-1)*8,(j-1)*8,:);
                mode = 5;
                [result_block,zz,cof] = intra_predictor(block,mode,qScale,prev_hori,prev_verti,0,diagonalvalue);
                result((i-1)*8+1:i*8,(j-1)*8+1:j*8,:) = result_block;
                zz_all = [zz_all,zz(:)'];
                pred_coeff = [pred_coeff,cof];

            else
                block = image((i-1)*8+1:i*8,(j-1)*8+1:j*8,:);
                prev_verti = result((i-1)*8,(j-1)*8+1:j*8,:);
                prev_verti2 = result((i-1)*8,j*8+1:(j+1)*8,:);
                prev_hori = result((i-1)*8+1:i*8,(j-1)*8,:);
                diagonalvalue = result((i-1)*8,(j-1)*8,:);
                mode = 4;
                [result_block,zz,cof] = intra_predictor(block,mode,qScale,prev_hori,prev_verti,prev_verti2,diagonalvalue);
                result((i-1)*8+1:i*8,(j-1)*8+1:j*8,:) = result_block;
                zz_all = [zz_all,zz(:)'];
                pred_coeff = [pred_coeff,cof];

            end
        end
    end   
recon_image = result;
dst = ZeroRunEnc_EoB(zz_all,1000);        
end




