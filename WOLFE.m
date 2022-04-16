function [satisfy, iWout] = WOLFE(x, al, d, f, g, c1, c2, iW)
    iWout = 0;
    satisfy = false;
    if f(x + al*d) <= f(x) + c1*g(x)'*d*al %satisfies WC1
        iWout = 1;
        if iW == 1 
            if g(x + al*d)'*d >= c2*g(x)'*d %satisfies WC2
                iWout = 2;
                satisfy = true;
            end
        elseif iW == 2 
            if abs(g(x + al*d)'*d) <= c2*abs(g(x)'*d)%satisfies SWC2
                iWout = 3;
                satisfy = true;
            end
        end
    end
end  