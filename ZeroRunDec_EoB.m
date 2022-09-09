function dst = ZeroRunDec_EoB(src,EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)
    % initial value
    idxDST = 1;
    i = 1;
    while 1

        if src(1,i) == EoB 
            % case 1: EOB
            if mod(idxDST,64) == 0
                residZeroOfBlock = 0;
            else
                residZeroOfBlock = 64 - mod(idxDST,64);
            end
                   
            dst(1,idxDST:idxDST+residZeroOfBlock) = 0;
            idxDST = idxDST+residZeroOfBlock+1;
            i = i+1;
        elseif src(1,i) ~= 0
            % case 2: nonzero integer
            dst(1,idxDST) = src(1,i);
            idxDST = idxDST + 1; 
            i = i+1;
        else
            % case 3: zero
            dst(1,idxDST:idxDST+src(1,i+1)) = src(1,i);
            idxDST = idxDST +src(1,i+1)+1;              
            i = i+2;
        end
        if i == numel(src)+1
            break
        end        
    end
end