function rgb = ictYCbCr2RGB(yuv)
% Input         : yuv (Original YCbCr image)
% Output        : rgb (RGB Image after transformation)
% YOUR CODE HERE
rgb(:,:,1) = 1*yuv(:,:,1) + 0*yuv(:,:,2) + 1.402*yuv(:,:,3);
rgb(:,:,2) = 1*yuv(:,:,1) -0.344*yuv(:,:,2) -0.714*yuv(:,:,3);
rgb(:,:,3) = 1*yuv(:,:,1) + 1.772*yuv(:,:,2) - 0*yuv(:,:,3);
end