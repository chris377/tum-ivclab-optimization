function rec_image = SSD_rec(ref_image, motion_vectors)
%  Input         : ref_image(Reference Image, YCbCr image)
%                  motion_vectors
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

[H_MV,W_MV,C]= size(motion_vectors);
for i = 1:H_MV
    for j= 1: W_MV
        [col,row] = ind2sub([9 9],motion_vectors(i,j));
        diff_col = col - 5;
        diff_row = row - 5;
        idy = (i-1)*8 + 1 + diff_row;
        idx = (j-1)*8 + 1 + diff_col;
      
        rec_block = ref_image(idy:idy+7,idx:idx+7,:);
        rec_image(i*8-7:i*8,j*8-7:j*8,:) = rec_block;
    end
end
end