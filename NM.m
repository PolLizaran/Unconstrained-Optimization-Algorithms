% Mètode de Newton amb longitud de pas = 1

function [xk, dk, alk, iWk, Hk] = NM(x, g, h, epsG, kmax) 
    xk = [x]; dk = []; alk = []; iWk =[];  Hk = [h(x)]; k = 1;
    al = 1;
    while norm(g(x)) > epsG & k < kmax
        p = -h(x)^-1 * g (x);
        x = x + al*p; 
        xk = [xk, x]; dk = [dk, p]; alk = [alk, 1]; iWk = [iWk, 4]; Hk(:,:,end+1) = [h(x)];
        k = k + 1;
    end
end