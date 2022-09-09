clc;clear all; close all; 
%% find out N
filep = mfilename('fullpath');
[pathstr,namestr]=fileparts(filep);
currentFolder = pathstr;
if filep(1,end-4) == '\'
    nextFolder = '\foreman20_40_RGB\';
else
    nextFolder = '/foreman20_40_RGB/';
end
fileDir = [currentFolder, nextFolder];
filePattern = [fileDir,'*.bmp'];
dirOutput = dir(filePattern); 
N = numel(dirOutput);
%% Initial value
rate = zeros(N,1);
PSNR = zeros(N,1);
%% For foreman0020.bmp image
lena_small = double(imread([fileDir,'lena_small.tif']));
lena_small_yuv = ictRGB2YCbCr(lena_small);
foreman0020 = double(imread([fileDir,'foreman0020.bmp']));
range = -1000:4000;

%% --------------------video optimized-------------------------------------------
fprintf('---------compute the mean bitrate and PSNR for video optimized-----------\n')
scales = [0.07 0.2 0.4 0.8 1.5 2 3 4 4.5]; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1

final_rate_optimized = zeros(numel(scales),1);
final_PSNR_optimized = zeros(numel(scales),1);
pracPixel = 2; % 1:voll pixel
               % 2:half pixel
               % 4:quarter pixel

for scaleIdx = 1 : numel(scales)

    rate = zeros(N,1); % reset the rate
    PSNR = zeros(N,1); % reset the PSNR  
    qScale   = scales(scaleIdx);

    %% use pmf of k_small to build and train huffman table
     
    [recon_image_small,k_small_videooptimized,pred_coeff_small_videoopmized]  = IntraEncode_new(lena_small_yuv, qScale);
    pmfRes    = stats_marg(k_small_videooptimized, range);   
    [ BinaryTree, HuffCode, BinCode, Codelengths] = buildHuffman(pmfRes );

    %% For foreman0020.bmp image
    foreman0020_yuv = ictRGB2YCbCr(foreman0020);
    [recon_image_videooptimized,k_videooptimized,pred_coeff_videoopmized] = IntraEncode_new(foreman0020_yuv, qScale);
    %% use trained table to encode k to get the bytestream
    bytestream  = enc_huffman_new( round(k_videooptimized(:))-range(1)+1, BinCode, Codelengths);
    bitPerPixel_foreman = (numel(bytestream)*8+numel(pred_coeff_videoopmized)) / (numel(foreman0020)/3);


    %% image reconstruction
    I_rec_rgb = ictYCbCr2RGB(recon_image_videooptimized);
    PSNR_foreman = calcPSNR(foreman0020, I_rec_rgb);

    % update rate(1) and PSNR(1) 
    rate(1) = bitPerPixel_foreman;
    PSNR(1) = PSNR_foreman;
    
    %% 2.frame to N.frame - Video Codec
    ref_im = I_rec_rgb;
    for frameIdx = 2:1:N
        if frameIdx == 2
            % train Huffmann code for MVs on first MV (between first and second frame)
            im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
            im_yuv = ictRGB2YCbCr(im);
            ref_im_yuv = ictRGB2YCbCr(ref_im);

            [err_im, MV] = MotionEstimationFractionalPixel(ref_im_yuv, im_yuv,pracPixel);
            dst = IntraEncode(err_im, qScale);  

            % build huffmann table using first residual error image
            pmfRes    = stats_marg(dst, range);    
            [ BinaryTreeRes, HuffCodeRes, BinCodeRes, CodelengthsRes] = buildHuffman(pmfRes );

            % build huffmann table using first motion vectors
            pmfMV    = stats_marg(MV, range);    
            [ BinaryTreeMV, HuffCodeMV, BinCodeMV, CodelengthsMV] = buildHuffman(pmfMV);            

        else
            im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
            im_yuv = ictRGB2YCbCr(im);
            ref_im_yuv = ictRGB2YCbCr(ref_im);
        
            [err_im, MV] = MotionEstimationFractionalPixel(ref_im_yuv, im_yuv,pracPixel);
            dst = IntraEncode(err_im, qScale);
        end

        % use trained table to encode k to get the bytestream
        bytestreamZeroRun  = enc_huffman_new( round(dst(:))-range(1)+1, BinCodeRes, CodelengthsRes);
        bytestreamMotionVec  = enc_huffman_new( round(MV(:))-range(1)+1, BinCodeMV, CodelengthsMV);
        rateZeroRun = (numel(bytestreamZeroRun)*8) / (numel(im)/3);     % rate1
        rateMotionVec = (numel(bytestreamMotionVec)*8) / (numel(im)/3); % rate2
        
        % image reconstruction
        zeroRun_rec  = double(reshape(dec_huffman_new (bytestreamZeroRun,  BinaryTreeRes, numel(dst)),  size(dst)))+range(1)-1;
        motionVec_rec  = double(reshape(dec_huffman_new (bytestreamMotionVec,  BinaryTreeMV, numel(MV)),  size(MV)))+range(1)-1;

        resi_rec = IntraDecode(zeroRun_rec, size(im) , qScale);
        recImg = SSD_rec_FractionalPixel(ref_im_yuv, motionVec_rec,pracPixel);
        loc_I_rec_ycbcr = recImg + resi_rec;
        I_rec_rgb = ictYCbCr2RGB(loc_I_rec_ycbcr);
    
        % calculate the rate and PSNR
        rate(frameIdx) = rateZeroRun + rateMotionVec; % Bit
        PSNR(frameIdx) = calcPSNR(im, I_rec_rgb); %PSNR
        ref_im = I_rec_rgb;

    
    end
    
    final_rate_optimized(scaleIdx) = mean(rate);
    final_PSNR_optimized(scaleIdx) = mean(PSNR);
     fprintf('QP Foreman(video codec): %.2f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, final_rate_optimized(scaleIdx),  final_PSNR_optimized(scaleIdx))

end
%% --------------------video original-------------------------------------------
fprintf('---------compute the mean bitrate and PSNR for video original-----------\n')
scales = [0.07 0.2 0.4 0.8 1.5 2 3 4 4.5]; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1

final_rate = zeros(numel(scales),1);
final_PSNR = zeros(numel(scales),1);

for scaleIdx = 1 : numel(scales)

    rate = zeros(N,1); % reset the rate
    PSNR = zeros(N,1); % reset the PSNR  
    qScale   = scales(scaleIdx);

    %% use pmf of k_small to build and train huffman table
    k_small  = IntraEncode(lena_small_yuv, qScale);
    pmfRes    = stats_marg(k_small, range);   
    [ BinaryTree, HuffCode, BinCode, Codelengths] = buildHuffman(pmfRes );

    %% For foreman0020.bmp image
    foreman0020_yuv = ictRGB2YCbCr(foreman0020);
    k        = IntraEncode(foreman0020_yuv, qScale);
    %% use trained table to encode k to get the bytestream
    bytestream  = enc_huffman_new( round(k(:))-range(1)+1, BinCode, Codelengths);
    bitPerPixel_foreman = (numel(bytestream)*8) / (numel(foreman0020)/3);


    %% image reconstruction
    k_rec  = double(reshape(dec_huffman_new (bytestream,  BinaryTree, numel(k)),  size(k)))+range(1)-1;
    I_rec_ycbcr = IntraDecode(k_rec, size(foreman0020),qScale);
    I_rec_rgb = ictYCbCr2RGB(I_rec_ycbcr);
    PSNR_foreman = calcPSNR(foreman0020, I_rec_rgb);

    % update rate(1) and PSNR(1) 
    rate(1) = bitPerPixel_foreman;
    PSNR(1) = PSNR_foreman;
    
    %% 2.frame to N.frame - Video Codec
    ref_im = I_rec_rgb;
    for frameIdx = 2:1:N
        if frameIdx == 2
            % train Huffmann code for MVs on first MV (between first and second frame)
            im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
            im_yuv = ictRGB2YCbCr(im);
            ref_im_yuv = ictRGB2YCbCr(ref_im);
        
            [err_im, MV] = MotionEstimation(ref_im_yuv, im_yuv);
    
            dst = IntraEncode(err_im, qScale);  

            % build huffmann table using first residual
            pmfRes    = stats_marg(dst, range);    
            [ BinaryTreeRes, HuffCodeRes, BinCodeRes, CodelengthsRes] = buildHuffman(pmfRes );

            % build huffmann table using first motion vectors
            pmfMV    = stats_marg(MV, range);    
            [ BinaryTreeMV, HuffCodeMV, BinCodeMV, CodelengthsMV] = buildHuffman(pmfMV);            

        else
            im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
            im_yuv = ictRGB2YCbCr(im);
            ref_im_yuv = ictRGB2YCbCr(ref_im);
        
            [err_im, MV] = MotionEstimation(ref_im_yuv, im_yuv);
    
            dst = IntraEncode(err_im, qScale);
        end

        % use trained table to encode k to get the bytestream
        bytestreamZeroRun  = enc_huffman_new( round(dst(:))-range(1)+1, BinCodeRes, CodelengthsRes);
        bytestreamMotionVec  = enc_huffman_new( round(MV(:))-range(1)+1, BinCodeMV, CodelengthsMV);
        rateZeroRun = (numel(bytestreamZeroRun)*8) / (numel(im)/3);     % rate1
        rateMotionVec = (numel(bytestreamMotionVec)*8) / (numel(im)/3); % rate2
        
        % image reconstruction
        zeroRun_rec  = double(reshape(dec_huffman_new (bytestreamZeroRun,  BinaryTreeRes, numel(dst)),  size(dst)))+range(1)-1;
        motionVec_rec  = double(reshape(dec_huffman_new (bytestreamMotionVec,  BinaryTreeMV, numel(MV)),  size(MV)))+range(1)-1;

        resi_rec = IntraDecode(zeroRun_rec, size(im) , qScale);
        recImg = SSD_rec(ref_im_yuv, motionVec_rec);
        loc_I_rec_ycbcr = recImg + resi_rec;
        I_rec_rgb = ictYCbCr2RGB(loc_I_rec_ycbcr);
    
        % calculate the rate and PSNR
        rate(frameIdx) = rateZeroRun + rateMotionVec; % Bit
        PSNR(frameIdx) = calcPSNR(im, I_rec_rgb); %PSNR
        ref_im = I_rec_rgb;

    
    end
    
    
    final_rate(scaleIdx) = mean(rate);
    final_PSNR(scaleIdx) = mean(PSNR);
     fprintf('QP Foreman(video codec): %.2f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, final_rate(scaleIdx),  final_PSNR(scaleIdx))

end


%% ---------------------still image optimized---------------------------------------------
fprintf('---------compute the mean bitrate and PSNR for still image optimized-----------\n')
scales = [0.15 0.3 0.7 1.0 1.5 3 5 7 10]; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1
final_rate_still_optimized = zeros(numel(scales),1);
final_PSNR_still_optimized = zeros(numel(scales),1);

for scaleIdx = 1 : numel(scales)
    rate = zeros(N,1); % reset the rate
    PSNR = zeros(N,1); % reset the PSNR
    qScale   = scales(scaleIdx);

    %% use pmf of k_small to build and train huffman table
    [recon_image_lena_small,k_small_optimized,pred_coeff_small]  = IntraEncode_new(lena_small_yuv, qScale);
    pmfRes_optimized    = stats_marg(k_small_optimized, range);   
    [ BinaryTree, HuffCode, BinCode, Codelengths] = buildHuffman(pmfRes_optimized );

    for frameIdx = 1:1:N
        im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
        im_yuv = ictRGB2YCbCr(im);
        [recon_image_foreman,k_optimized,pred_coeff_foreman] = IntraEncode_new(im_yuv, qScale);
        %% use trained table to encode k to get the bytestream
        bytestream  = enc_huffman_new( round(k_optimized(:))-range(1)+1, BinCode, Codelengths);
        bitPerPixel_foreman = (numel(bytestream)*8+numel(pred_coeff_foreman)) / (numel(im_yuv)/3);
    
        %% image reconstruction
        I_rec_rgb = ictYCbCr2RGB(recon_image_foreman);
        PSNR_foreman = calcPSNR(im, I_rec_rgb);
    
        % update rate(1) and PSNR(1) 
        rate(frameIdx) = bitPerPixel_foreman; % Bit
        PSNR(frameIdx) = PSNR_foreman; %PSNR
    
    end
    
    
    final_rate_still_optimized(scaleIdx) = mean(rate);
    final_PSNR_still_optimized(scaleIdx) = mean(PSNR);
    fprintf('QP Foreman(Still image codec): %.2f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, final_rate_still_optimized(scaleIdx),  final_PSNR_still_optimized(scaleIdx))
end



%% ---------------------still image origin---------------------------------------------
fprintf('---------compute the mean bitrate and PSNR for still image original-----------\n')
scales = [0.15 0.3 0.7 1.0 1.5 3 5 7 10]; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1
final_rate_still = zeros(numel(scales),1);
final_PSNR_still = zeros(numel(scales),1);

for scaleIdx = 1 : numel(scales)
    rate = zeros(N,1); % reset the rate
    PSNR = zeros(N,1); % reset the PSNR
    qScale   = scales(scaleIdx);

    %% use pmf of k_small to build and train huffman table
    k_small  = IntraEncode(lena_small_yuv, qScale);
    pmfRes    = stats_marg(k_small, range);   
    [ BinaryTree, HuffCode, BinCode, Codelengths] = buildHuffman(pmfRes );

    for frameIdx = 1:1:N
        im = double(imread([fileDir,dirOutput(frameIdx).name])) ;
        im_yuv = ictRGB2YCbCr(im);
        k        = IntraEncode(im_yuv, qScale);
        %% use trained table to encode k to get the bytestream
        bytestream  = enc_huffman_new( round(k(:))-range(1)+1, BinCode, Codelengths);
        bitPerPixel_foreman = (numel(bytestream)*8) / (numel(im_yuv)/3);
    
        %% image reconstruction
        k_rec  = double(reshape(dec_huffman_new (bytestream,  BinaryTree, numel(k)),  size(k)))+range(1)-1;
        loc_I_rec_ycbcr= IntraDecode(k_rec, size(im_yuv),qScale);
        I_rec_rgb = ictYCbCr2RGB(loc_I_rec_ycbcr);
        PSNR_foreman = calcPSNR(im, I_rec_rgb);
    
        % update rate(1) and PSNR(1) 
        rate(frameIdx) = bitPerPixel_foreman; % Bit
        PSNR(frameIdx) = PSNR_foreman; %PSNR
    
    end
    
    
    final_rate_still(scaleIdx) = mean(rate);
    final_PSNR_still(scaleIdx) = mean(PSNR);
    fprintf('QP Foreman(Still image codec): %.2f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, final_rate_still(scaleIdx),  final_PSNR_still(scaleIdx))
end








%% Plot 
figure(1);
title('D-R curve')
xlabel('bpp') 
ylabel('PSNR [dB]')

hold on
plot(final_rate_optimized, final_PSNR_optimized,'--*')
plot(final_rate, final_PSNR,'--*')
plot(final_rate_still_optimized, final_PSNR_still_optimized,'--*')
plot(final_rate_still, final_PSNR_still,'--*')
legend('video_optimized','video_original','image_optimized','image_original')
hold off
