function zze = ZeroRunEnc_EoB(zz, EoB)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)
    % initial value
    idxZZE = 1;
    zeroCounter = -1;
    for i = 1: numel(zz)
        if zz(1,i) ~= 0
            if zeroCounter ~= -1
                zze(1,idxZZE) = zeroCounter;
                idxZZE = idxZZE+1;
                zeroCounter=-1;
                
                zze(1,idxZZE) = zz(1,i);
                idxZZE = idxZZE+1;                
            else
                zze(1,idxZZE) = zz(1,i);
                idxZZE = idxZZE+1;                           
            end
        else
            if zz(1,i) == 0 && zeroCounter == -1
                zze(1,idxZZE) = 0;
                idxZZE = idxZZE+1;
                zeroCounter = zeroCounter+1;
            else
                zeroCounter = zeroCounter+1;
            end
            
            % If the last number of the whole matrix or the last number of 8*8 block is 0, then the encoder should end up with EOB.
            if i == numel(zz) || ((i/64) == round(i/64)) 
                zeroCounter = -1;
                zze(1,idxZZE-1) = EoB;
            end
        end
    end
end
