function gamma_k = gamma_(k,R_pos)
    gamma_k=0;
    delta_moy= mean(diff(R_pos));
    for n=1:length(R_pos)-k
        gamma_k=gamma_k+ ((R_pos(n+k+1)-R_pos(n+k)-delta_moy)*(R_pos(n+1)-R_pos(n)))/(length(R_pos)-k-1);
    end
end
