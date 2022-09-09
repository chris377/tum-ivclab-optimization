function [error,predict] = predict_mode2(image_block,prev_hori,prev_verti,diagonalvalue)
    predict = zeros(size(image_block));
    [x,y,z] = size(image_block);
    for k = 1:z
        predict(:,:,k) = diag(diagonalvalue(:,:,k)*ones(1,x));
    end
    
    for k = 1:z
        for i = 1:x
            predict(i,i+1:end,k) = prev_verti(1,1:end-i,k);
        end
        for i = 1:y
            predict(i+1:end,i,k) = prev_hori(1:end-i,1,k);
        end
    end
    error = image_block - predict;

end