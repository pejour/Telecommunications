function [bits_reconstruit] = demapping_4_ASK(decision)
    n = length(decision);
    bits_reconstruit = zeros(1, 2*n);
    for j = 1:n
       if (decision(j) == -3)
           bits_reconstruit(2*j) = 0;
           bits_reconstruit(2*j - 1) = 0;
       elseif (decision(j) == -1)
           bits_reconstruit(2*j) = 1;
           bits_reconstruit(2*j - 1) = 0;
       elseif (decision(j) == 1)
           bits_reconstruit(2*j) = 1;
           bits_reconstruit(2*j - 1) = 1;
       else
           bits_reconstruit(2*j) = 0;
           bits_reconstruit(2*j - 1) = 1;
       end
    end
end

