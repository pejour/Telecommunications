function [decision] = decision_8_PSK(ech)
    n = length(ech);
    for j = 1:n
        if (0 <= angle(ech(j)) && angle(ech(j)) < (pi / 4))
            decision(j) = exp(1i * (pi/8));
        elseif ((pi/4) <= angle(ech(j)) && angle(ech(j)) < (pi / 2))
            decision(j) = exp(1i * (((2 * pi )/ 8) + (pi / 8)));
        elseif ((pi / 2) <= angle(ech(j)) && angle(ech(j)) < (3 * pi / 4))
            decision(j) = exp(1i * (((2 * pi )/ 8) * 2 + (pi / 8)));
        elseif ((3 * pi /4) <= angle(ech(j)) && angle(ech(j)) < pi)
            decision(j) = exp(1i * (((2 * pi )/ 8) * 3 + (pi / 8)));
        elseif (-pi <= angle(ech(j)) && angle(ech(j)) < (-3 * pi / 4))
            decision(j) = exp(1i * (((2 * pi )/ 8) * 4 + (pi / 8)));
        elseif ((-3 * pi / 4) <= angle(ech(j)) && angle(ech(j)) < (-pi / 2))
            decision(j) = exp(1i * (((2 * pi )/ 8) * 5 + (pi / 8)));
        elseif ((-pi / 2) <= angle(ech(j)) && angle(ech(j)) < (-pi / 4))
            decision(j) = exp(1i * (((2 * pi )/ 8) * 6 + (pi / 8)));
        else
            decision(j) = exp(1i * (((2 * pi )/ 8) * 7 + (pi / 8)));
        end
    end
end

