% Funció que usa el Mètode del Gradient: d = -g(x). La direcció de moviment
% es pot calcular amb un ELS o BLS

function [xk, dk, alk, iWk] = GM(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW)

    xk = [x]; dk = []; alk = []; iWk =[]; k = 1;
    while norm(g(x)) >= epsG & k < kmax
        d = -g(x); dk = [dk, d];
        [al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW);
        x = x + d*al;
        xk = [xk, x];  alk = [alk, al]; iWk = [iWk, iWout]; 
        k = k + 1;
    end
    disp(k);
end