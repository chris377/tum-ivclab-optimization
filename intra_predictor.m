function  [result,zz,cof] = intra_predictor(block,mode,qScale,prev_hori,prev_verti,prev_verti2,diagvalue)
    switch mode
        case 1
            coeff = DCT8x8(block);
            quant = Quant8x8(coeff, qScale);
            zz    = ZigZag8x8(quant);
            recon_block = DeQuant8x8(quant,qScale);
            result = IDCT8x8(recon_block);
            cof = "";
            
        case 2
            [error,predict] = predict_mode0(block,prev_hori,prev_verti);
            coeff = DCT8x8(error);
            quant = Quant8x8(coeff,qScale);  
            zz    = ZigZag8x8(quant);
            recon_block = DeQuant8x8(quant,qScale);
            result = IDCT8x8(recon_block)+predict;
            cof = '00';
        case 3
            [error1,predict1] = predict_mode1(block,prev_hori,prev_verti);
            [error2,predict2] = predict_mode4(block,prev_verti,prev_verti2);
            [error,predict,cof] = cal_sad3(block,error1,predict1,error2,predict2);
            coeff = DCT8x8(error);
            quant = Quant8x8(coeff,qScale);  
            zz    = ZigZag8x8(quant);
            recon_block = DeQuant8x8(quant,qScale);
            result = IDCT8x8(recon_block)+predict;
            %cof = '01';
        case 4
            [error1,predict1] = predict_mode0(block,prev_hori,prev_verti);
            [error2,predict2] = predict_mode1(block,prev_hori,prev_verti);
            [error3,predict3] = predict_mode2(block,prev_hori,prev_verti,diagvalue);
            [error4,predict4] = predict_mode3(block,prev_hori,prev_verti);
            [error5,predict5] = predict_mode4(block,prev_verti,prev_verti2);
            [error,predict,cof] = cal_sad1(block,error1,predict1,error2,predict2,error3,predict3,error4,predict4,error5,predict5);
            coeff = DCT8x8(error);
            quant = Quant8x8(coeff,qScale);  
            zz    = ZigZag8x8(quant);
            recon_block = DeQuant8x8(quant,qScale);
            result = IDCT8x8(recon_block)+predict;

        case 5
            [error1,predict1] = predict_mode0(block,prev_hori,prev_verti);
            [error2,predict2] = predict_mode1(block,prev_hori,prev_verti);
            [error3,predict3] = predict_mode2(block,prev_hori,prev_verti,diagvalue);
            [error4,predict4] = predict_mode3(block,prev_hori,prev_verti);
            [error,predict,cof] = cal_sad2(block,error1,predict1,error2,predict2,error3,predict3,error4,predict4);
            coeff = DCT8x8(error);
            quant = Quant8x8(coeff,qScale);  
            zz    = ZigZag8x8(quant);
            recon_block = DeQuant8x8(quant,qScale);
            result = IDCT8x8(recon_block)+predict;
            
    end
    result(result(:,:,1) > 255) = 255;
    result(result(:,:,1) < 0  ) = 0;
end