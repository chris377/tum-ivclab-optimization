function [error,predict] = predict_mode0(image_block,prev_hori,prev_verti)
    [x,y,z] = size(image_block);
    predict = zeros(size(image_block));
    for i = 1:z
        predict(:,:,i) = repmat(prev_hori(:,:,i),1,y);
    end
    error = image_block - predict;

end
