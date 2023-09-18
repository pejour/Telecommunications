function [bits] = demapping_8_PSK(decision)
    n = length(decision);
    bits = zeros(1,3 * n);
    for j = 1:n
        if (decision(j) == exp(1i * (pi/8)))   % On reconstruit le bit 111
            bits(3*j) = 1;
            bits(3*j - 1) = 1;
            bits(3*j - 2) = 1;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) + (pi / 8)))) % bit 011
            bits(3*j) = 0;
            bits(3*j - 1) = 1;
            bits(3*j - 2) = 1;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) * 2 + (pi / 8)))) % bit 010
            bits(3*j) = 0;
            bits(3*j - 1) = 1;
            bits(3*j - 2) = 0;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) * 3 + (pi / 8))))  % bit 110
            bits(3*j) = 1;
            bits(3*j - 1) = 1;
            bits(3*j - 2) = 0;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) * 4 + (pi / 8))))  % bit 100
            bits(3*j) = 1;
            bits(3*j - 1) = 0;
            bits(3*j - 2) = 0;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) * 5 + (pi / 8))))  % bit 000
            bits(3*j) = 0;
            bits(3*j - 1) = 0;
            bits(3*j - 2) = 0;
        elseif (decision(j) == exp(1i * (((2 * pi )/ 8) * 6 + (pi / 8))))  % bit 001
            bits(3*j) = 0;
            bits(3*j - 1) = 0;
            bits(3*j - 2) = 1;
        else                        % bit 101
            bits(3*j) = 1;
            bits(3*j - 1) = 0;
            bits(3*j - 2) = 1;
        end
    end
end

