  
ext=[];
RR= ecg_g2(position(1):position(1+1));
    for j=position(1):round(0.7*(position(1+1)-position(1)))+position(1)
        if ecg_g2(j) >0 && ecg_g2(j+1)<0 % tu prend j 
           ext=[ext j]; % Tous les maximum locaux
        end
    end