
function interpolated_image = QuartelPixelInterpolation(image)
%  Input         : image(YCbCr image with size H*W*C)
%
%  Output        : interpolated_image(YCbCr image with size
%                  (4*(H-1)+1)*(4*(W-1)+1)*C

    % Initial value
    [H, W, C] = size(image);
    halfInterplated_image = zeros(2*(H-1)+1,2*(W-1)+1,C);
    interpolated_image = zeros(4*(H-1)+1,4*(W-1)+1,C);


    if 0 % with round
        % 6-tap-filtering
        % first step: create a 1/2 piexl interpolaton 
        % Example along the rows (colums same)
        % AaBbCcDdEeFF
        % c = round((A-5*B+20*C+20*D-5*E+F)/32)        
        halfInterplated_image(1:2:end,1:2:end,:) = image;
        halfInterplated_image((2*3):2:end-(2*3-1),:,:) = round(( ...
            halfInterplated_image(1:2:end-10,:,:)-5*halfInterplated_image(3:2:end-8,:,:)...
            +20*halfInterplated_image(5:2:end-6,:,:)+20*halfInterplated_image(7:2:end-4,:,:)...
            -5*halfInterplated_image(9:2:end-2,:,:)+halfInterplated_image(11:2:end,:,:)...
            )/32); % interplation along the rows
        halfInterplated_image(:,(2*3):2:end-(2*3-1),:) = round(( ...
            halfInterplated_image(:,1:2:end-10,:)-5*halfInterplated_image(:,3:2:end-8,:)...
            +20*halfInterplated_image(:,5:2:end-6,:)+20*halfInterplated_image(:,7:2:end-4,:)...
            -5*halfInterplated_image(:,9:2:end-2,:)+halfInterplated_image(:,11:2:end,:)...
            )/32); % interpolation along the columns
    
        halfInterplated_image(:,[2,4,end-3,end-1],:) = round(( ...
            halfInterplated_image(:,[2,4,end-3,end-1]-1,:)+...
            halfInterplated_image(:,[2,4,end-3,end-1]+1,:)...
            )/2); % interplation along the 1 2 end-1 end rows
        halfInterplated_image([2,4,end-3,end-1],:,:) = round(( ...
            halfInterplated_image([2,4,end-3,end-1]-1,:,:)+...
            halfInterplated_image([2,4,end-3,end-1]+1,:,:)...
            )/2); % interplation along the 1 2 end-1 end columns
    
        % 6-tap-filtering
        % second step: create a 1/4 piexl interpolaton by interpolation
        % of the adjacent 1/2 pixel
        interpolated_image(1:2:end,1:2:end,:) = halfInterplated_image;
        interpolated_image(1:2:end,2:2:end-1,:) = round((...
            interpolated_image(1:2:end,1:2:end-2,:)+...
            interpolated_image(1:2:end,3:2:end,:))/2); % along row
        interpolated_image(2:2:end-1,1:2:end,:) = round((...
            interpolated_image(1:2:end-2,1:2:end,:)+...
            interpolated_image(3:2:end,1:2:end,:))/2); % along columns
        interpolated_image(2:2:end-1,2:4:end-3,:) = round((...
            interpolated_image(3:2:end,1:4:end-4,:)+...
            interpolated_image(1:2:end-2,3:4:end-2,:))/2); % along diagonal
        interpolated_image(2:2:end-1,4:4:end-1,:) = round((...
            interpolated_image(3:2:end,5:4:end,:)+...
            interpolated_image(1:2:end-2,3:4:end-2,:))/2); % along diagonal 

    else % without round
        % 6-tap-filtering
        % first step: create a 1/2 piexl interpolaton 
        % Example along the rows (colums same)
        % AaBbCcDdEeFF
        % c = round((A-5*B+20*C+20*D-5*E+F)/32)            
        halfInterplated_image(1:2:end,1:2:end,:) = image;
        halfInterplated_image((2*3):2:end-(2*3-1),:,:) = ( ...
            halfInterplated_image(1:2:end-10,:,:)-5*halfInterplated_image(3:2:end-8,:,:)...
            +20*halfInterplated_image(5:2:end-6,:,:)+20*halfInterplated_image(7:2:end-4,:,:)...
            -5*halfInterplated_image(9:2:end-2,:,:)+halfInterplated_image(11:2:end,:,:)...
            )/32; % interplation along the rows
        halfInterplated_image(:,(2*3):2:end-(2*3-1),:) = ( ...
            halfInterplated_image(:,1:2:end-10,:)-5*halfInterplated_image(:,3:2:end-8,:)...
            +20*halfInterplated_image(:,5:2:end-6,:)+20*halfInterplated_image(:,7:2:end-4,:)...
            -5*halfInterplated_image(:,9:2:end-2,:)+halfInterplated_image(:,11:2:end,:)...
            )/32; % interpolation along the columns
    
        halfInterplated_image(:,[2,4,end-3,end-1],:) = ( ...
            halfInterplated_image(:,[2,4,end-3,end-1]-1,:)+...
            halfInterplated_image(:,[2,4,end-3,end-1]+1,:)...
            )/2; % interplation along the 1 2 end-1 end rows
        halfInterplated_image([2,4,end-3,end-1],:,:) = ( ...
            halfInterplated_image([2,4,end-3,end-1]-1,:,:)+...
            halfInterplated_image([2,4,end-3,end-1]+1,:,:)...
            )/2; % interplation along the 1 2 end-1 end columns
    
        % 6-tap-filtering
        % first step: create a 1/4 piexl interpolaton by interpolation
        % of the adjacent 1/2 pixel
        interpolated_image(1:2:end,1:2:end,:) = halfInterplated_image;
        interpolated_image(1:2:end,2:2:end-1,:) = (...
            interpolated_image(1:2:end,1:2:end-2,:)+...
            interpolated_image(1:2:end,3:2:end,:))/2; % along row
        interpolated_image(2:2:end-1,1:2:end,:) = (...
            interpolated_image(1:2:end-2,1:2:end,:)+...
            interpolated_image(3:2:end,1:2:end,:))/2; % along columns
        interpolated_image(2:2:end-1,2:4:end-3,:) = (...
            interpolated_image(3:2:end,1:4:end-4,:)+...
            interpolated_image(1:2:end-2,3:4:end-2,:))/2; % along diagonal
        interpolated_image(2:2:end-1,4:4:end-1,:) = (...
            interpolated_image(3:2:end,5:4:end,:)+...
            interpolated_image(1:2:end-2,3:4:end-2,:))/2; % along diagonal        
    end    
 
end