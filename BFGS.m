% Calcula el punt òptim mitjançant el mètode BFGS d'aproximació de la
% matriu Hessiana. 

function [xk, dk, alk, iWk, betak, Hk] = BFGS(x, f, g, h, epsG, kmax, almax, ...
                                              almin, rho, c1, c2, iW)

    n = size(x); n1 = n(1); I = eye(n1); H = I; % Inicialització de la matriu H
    xk = [x]; dk = []; alk = []; iWk =[]; betak = []; Hk = [H]; k = 1; 
    while norm(g(x)) > epsG & k < kmax
        d = -H * g(x);
        [al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW);
        x = x + d*al; 
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al]; iWk = [iWk, iWout];

        s = xk(:, end) - xk(:, end - 1); % Equació secant
        y = g(xk(:,end)) - g(xk(:,end - 1)); % Diferència entre gradients
        p = 1/(y'*s); 
        H = (I - p*s*y')*H*(I - p*y*s') + p*s*s'; Hk(:,:,end+1) = [H];
        k = k + 1; 
    end
end
