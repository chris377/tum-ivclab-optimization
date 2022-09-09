function [error,predict] = predict_mode4(image_block,prev_verti1,prev_verti2)
    [x,y,z] = size(image_block);
    predict = zeros(size(image_block));

    for k = 1:z
        for i = 1:x-1
            predict(i,1:x-i,k) = prev_verti1(1,i+1:end,k);
        end
        for j = 1:y
            predict(y-j+1:end,j,k) = prev_verti2(1,1:j,k);
        end
    end
    error = image_block - predict;

end