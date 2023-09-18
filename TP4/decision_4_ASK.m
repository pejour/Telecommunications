function [decision] = decision_4_ASK(ech)
    n = length(ech);
    decision = zeros(1, n);
    for i = 1:n
        if (ech(i) <= -2)
            decision(i) = -3;
        elseif (-2 < ech(i) && ech(i) <= 0)
            decision(i) = -1;
        elseif (0 < ech(i) && ech(i) <= 2)
            decision(i) = 1;
        else
            decision(i) = 3;
        end
    end
end

