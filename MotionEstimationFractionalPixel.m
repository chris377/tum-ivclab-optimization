%--------------------------------------------------------------------
function [err_im, MV] = MotionEstimationFractionalPixel(ref_image, image,N)
%  Input         : ref_image(Reference Image, size: height x width x channel)
%                  image (Current Image, size: height x width)
%                  N(fractional-pel accuracy 1/N)
%  Output        : MV (Motion Vector Indices, size: (height/8) x (width/8) x 1 )
%                  err_im(the difference of each pixel, size: height x width x channel)
[H,W,C] = size(image);
err_im = zeros(H,W,C);
recImg = zeros(H,W,C);
MV = zeros(H/8,W/8);

% To get image block of fractional pixel
switch N
    case 1
        img_block_fracPixel = ref_image;
    case 2
        % USE 1/2 interpolation based on 6-tap filtering
        halfInterplated_image = HalfPixel6TapInterpolation(ref_image);
        img_block_fracPixel = halfInterplated_image;        
    case 4

        % USE 1/4 interpolation based on 6-tap filtering
        quartelInterplated_image = QuartelPixelInterpolation(ref_image);
        img_block_fracPixel = quartelInterplated_image;   
end


% estimate motion vector
for i = 1:8: H
    for j = 1:8:W
        img_block = image(i:i+7,j:j+7,1); % luminance component
        
        idx_ref = 1;
        min_mse_block = inf;

        min_mse_idx = 1;
        for i_ref = N*(i-1)+1-4*N:1:N*(i-1)+1+4*N
            for j_ref = N*(j-1)+1-4*N:1:N*(j-1)+1+4*N
                if i_ref>0 && j_ref>0 && i_ref<=(N*(H-1)+1-(N*8-1)) &&  j_ref<=(N*(W-1)+1-(N*8-1))
                    if img_block_fracPixel(i_ref,j_ref)>=0
                        ref_block = img_block_fracPixel(i_ref:N:i_ref+(N*8-1),j_ref:N:j_ref+(N*8-1),1); % luminance components
    
    
                        diff_block = img_block - ref_block;  % luminance components
                        error_block = diff_block.^2;
                        mse_block = sum(error_block(:));
                        recImgBlock = img_block_fracPixel(i_ref:N:i_ref+(N*8-1),j_ref:N:j_ref+(N*8-1),:);
                        if mse_block <  min_mse_block 
                            min_mse_block = mse_block;
                            min_mse_idx = idx_ref;
                          min_recImgBlock = recImgBlock;
                        end 
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

