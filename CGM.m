% Mètode del Gradient Conjugat que calcula les betes mitjançant FR o PR

function [xk, dk, alk, iWk, betak] = CGM(x, f, g, h, epsG, kmax, almax, almin, rho, c1, ...
                                         c2, iW, icg, irc, nu)

    xk = [x]; d = -g(x); dk = [d]; alk = []; iWk =[]; betak = [0]; k = 1;
    d_positiva = false;
    while norm(g(x)) > epsG & k < kmax
        %[al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW);
        [al, iWout] = BLS_DC(x, d, f, g, h, almax, almin, rho, c1, c2, iW);
        x = x + d*al; 
        beta = BETA_CHOICE(x, xk, g, icg, irc, nu, k); % icg = 1: Fletcher-Reeves, icg = 2: Polak-Ribière
        d = -g(x) + beta*dk(:,end);
        if (g(x)'*d > 0) 
            d_positiva = true;
        end
        dk = [dk, d]; alk = [alk, al]; iWk = [iWk, iWout]; betak = [betak, beta]; 
        xk = [xk, x]; %no es pot afegir abans de trobar la beta, sinó estaríem buscant el gradient sobre k-1, que seria x
        k = k + 1;
    end
    if d_positiva == true
        fprintf("\n \n !!!HI HA DIRECCIONS D'ASCENS POSITIVES!!! \n");
    end
end
