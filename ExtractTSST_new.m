function Ex = ExtractTSST_new( Tx, IdealGD, fs,singleside,direction)
    [L,N] = size(Tx);
    Temp = zeros(L,N);
    dt = 1/fs;
    if direction==1 %1代表时间负反向
      for ptr = 1:L
        upboard = min( round(IdealGD(ptr)/dt) + singleside+1,N);
        downboard = 1;
        Temp(ptr,downboard:upboard) = 1;
      end
 
    elseif direction==2 %2代表时间正方向
      for ptr = 1:L
        upboard = N;
        downboard = max( round(IdealGD(ptr)/dt) - singleside+1,1);
        Temp(ptr,downboard:upboard) = 1;
     end
    else %3 为带通  
        for ptr = 1:L
        upboard = min( round(IdealGD(ptr)/dt) + singleside+1,N);
        downboard = max( round(IdealGD(ptr)/dt) - singleside+1,1);
        Temp(ptr,downboard:upboard) = 1;
       end
    end
      Ex = Tx .* Temp;  
end