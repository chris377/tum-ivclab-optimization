function [error,predict,cof] = cal_sad3(image_block,error1,predict1,error2,predict2)
    
    if sum(abs(image_block-predict1)) < sum(abs(image_block-predict2))
        error = error1;
        predict = predict1;
        cof = "00";
    else
        error = error2;
        predict = predict2;
        cof = "01";
    end

end
    