function [symboles] = mapping_8_PSK(bits)
    n = length(bits);
    symboles = zeros(1, floor(n/3));
    % on génère les symboles dans le plan complexe, en s'aidant du cercle
    % trigonométrique
    k = 1;
    for j = (1:3:n)
        if (bits(j) == 1  && bits(j+1) == 1 && bits(j+2) == 1)    % bit 111
            symboles(k) = exp(1i * (pi/8));
            k = k + 1;
        elseif (bits(j) == 1  && bits(j+1) == 1 && bits(j+2) == 0) % bit 110
            symboles(k) = exp(1i * (((2 * pi )/ 8) + (pi / 8)));
            k = k + 1;
        elseif (bits(j) == 0  && bits(j+1) == 1 && bits(j+2) == 0)  % bit 010
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 2 + (pi / 8)));
            k = k + 1;
        elseif (bits(j) == 0  && bits(j+1) == 1 && bits(j+2) == 1)  % bit 011
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 3 + (pi / 8)));
            k = k + 1;
        elseif (bits(j) == 0  && bits(j+1) == 0 && bits(j+2) == 1)  % bit 001
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 4 + (pi / 8)));
            k = k + 1;
        elseif (bits(j) == 0  && bits(j+1) == 0 && bits(j+2) == 0)  % bit 000
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 5 + (pi / 8)));
            k = k + 1;
        elseif (bits(j) == 1  && bits(j+1) == 0 && bits(j+2) == 0)  % bit 100
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 6 + (pi / 8)));
            k = k + 1;
        else
            % dernier symbole non généré correspondant au bit 101
            symboles(k) = exp(1i * (((2 * pi )/ 8) * 7 + (pi / 8))); 
            k = k + 1;
        end
    end
end        
    

