function [error,predict] = predict_mode3(image_block,prev_hori,prev_verti)
    predict = ones(8,8,3);
    predict(:,:,1) = mean([prev_hori(:,:,1);prev_verti(:,:,1)']);
    predict(:,:,2) = mean([prev_hori(:,:,2);prev_verti(:,:,2)']);
    predict(:,:,3) = mean([prev_hori(:,:,3);prev_verti(:,:,3)']);
    error = image_block - predict;
end