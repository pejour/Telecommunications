function [symboles] = mapping_4_ASK(bits)
    n = length(bits);
    symboles = zeros(1, fix(n/2));
    for i = 1:fix(n/2)
        if (bits(2*i) == 0 && bits(2*i - 1) == 0)
            symboles(i) = -3;
        elseif (bits(2*i) == 1 && bits(2*i - 1) == 0)
            symboles(i) = -1;
        elseif (bits(2*i) == 1 && bits(2*i - 1) == 1)
            symboles(i) = 1;
        else
            symboles(i) = 3;
        end
    end
end

