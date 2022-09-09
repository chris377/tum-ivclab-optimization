function MSE = calcMSE(Image, recImage)
    W = size(Image,1);  % The width of the image
    H = size(Image,2);  % The height of the image
    C = size(Image,3);  % Number of the channels of the image
    
    SumSquaredError = 0; % intial value
    sumTest = 0;
    for i = 1:W
        for j = 1:H
            for k = 1:C
                SumSquaredError = SumSquaredError + (double(Image(i,j,k))-double(recImage(i,j,k)))^2;
                sumTest = sumTest+1;
                
            end
        end
    end
    
    MSE = double(SumSquaredError)/(W*H*C);
end