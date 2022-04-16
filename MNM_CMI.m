% Mètode que usa aproximació de la Hessiana a partir de la factorització de
% Cholesky. 

function [xk, dk, alk, iWk, Hk, tauk] = ...
          MNM_CMI(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW)
    xk = [x]; dk = []; alk = []; iWk =[];  Hk = [h(x)]; tauk = []; k = 1;
    while norm(g(x)) > epsG & k < kmax
        [B, tau] = B_MODIF_CMI(x, h);
        d = -B^-1 * g(x);
        [al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW); 
        x = x + al*d; 
        Hk(:,:,end+1) = [B];
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al]; iWk = [iWk, iWout]; tauk = [tauk, tau];
        k = k + 1;
    end
end