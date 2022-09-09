 function [err_im, MV] = MotionEstimation(ref_image, image)
%  Input         : ref_image(Reference Image, size: height x width x channel)
%                  image (Current Image, size: height x width)
% 
%  Output        : MV (Motion Vector Indices, size: (height/8) x (width/8) x 1 )
%                  err_im(the difference of each pixel, size: height x width x channel)
[H,W,C] = size(image);
err_im = zeros(H,W,C);
recImg = zeros(H,W,C);
MV = zeros(H/8,W/8);

for i = 1:8: H
    for j = 1:8:W
        img_block = image(i:i+7,j:j+7,1); % luminance component
        
        idx_ref = 1;
        min_mse_block = inf;

        min_mse_idx = 1;
        for i_ref = i-4:1:i+4
            for j_ref = j-4:1:j+4
                if i_ref>0 && j_ref>0 && i_ref<=(H-7) &&  j_ref<=(W-7)
                    ref_block = ref_image(i_ref:i_ref+7,j_ref:j_ref+7,1); % luminance components
                    diff_block = img_block - ref_block;  % luminance components
                    error_block = diff_block.^2;
                    mse_block = sum(error_block(:));
                   recImgBlock = ref_image(i_ref:i_ref+7,j_ref:j_ref+7,:);
                    if mse_block <  min_mse_block 
                        min_mse_block = mse_block;
                        min_mse_idx = idx_ref;
                      min_recImgBlock = recImgBlock;
                    end 
                end
                idx_ref= idx_ref+1;
            end
        end
        MV(floor(i/8)+1,floor(j/8)+1) = min_mse_idx;
        recImg(i:i+7,j:j+7,:) = min_recImgBlock;
    end
end
err_im = image - recImg;

end