function [error,predict,cof] = cal_sad1(image_block,error1,predict1,error2,predict2,error3,predict3,error4,predict4,error5,predict5)
    
    if sum(abs(image_block-predict1)) < sum(abs(image_block-predict2))
        error = error1;
        predict = predict1;
        cof = "00";
    else
        error = error2;
        predict = predict2;
        cof = "01";
    end

    if sum(abs(image_block-predict3))< sum(abs(image_block-predict))
        error = error3;
        predict = predict3;
        cof = "10";
    end

    if sum(abs(image_block-predict4))< sum(abs(image_block-predict))
        error = error4;
        predict = predict4;
        cof = "110";
    end

    if sum(abs(image_block-predict5))< sum(abs(image_block-predict))
        error = error5;
        predict = predict5;
        cof = "111";
    end

end
    