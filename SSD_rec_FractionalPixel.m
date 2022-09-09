%---------------------------------------------------------------------
function rec_image = SSD_rec_FractionalPixel(image, motion_vectors,N)
%  Input         : image(Reference Image, YCbCr image)
%                  motion_vectors
%                  N(fractional-pel accuracy 1/N)
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

[H_MV,W_MV,C]= size(motion_vectors);
[H,W,C] = size(image);


% To get image block of halb pixel 
switch N
    case 1
        img_block_fracPixel = image;
    case 2
        % USE 1/2 interpolation based on 6-tap filtering
        halfInterplated_image = HalfPixel6TapInterpolation(image);
        img_block_fracPixel = halfInterplated_image;          
    case 4
        % USE 1/4 interpolation based on 6-tap filtering
        quartelInterplated_image = QuartelPixelInterpolation(image);
        img_block_fracPixel = quartelInterplated_image;           
end

% estimate motion vector
for i = 1:H_MV
    for j= 1: W_MV
        [col,row] = ind2sub([N*8+1 N*8+1],motion_vectors(i,j));
        diff_col = col - (N*4+1);
        diff_row = row - (N*4+1);
        idy = (i-1)*N*8 + 1 + diff_row;
        idx = (j-1)*N*8 + 1 + diff_col;

        rec_block = img_block_fracPixel(idy:N:idy+(N*8-1),idx:N:idx+(N*8-1),:);
        rec_image(i*8-7:i*8,j*8-7:j*8,:) = rec_block;
    end
end
end