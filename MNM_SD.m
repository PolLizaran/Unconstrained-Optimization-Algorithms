% Mètode de Newton Modificat que realiza una variació de la matriu Hessiana
% usant la descomposició en valors singulars

function [xk, dk, alk, iWk, Hk] = ...
          MNM_SD(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW, delta)
    vaps_negatius = false;
    xk = [x]; dk = []; alk = []; iWk =[]; Hk = []; k = 1;
    while norm(g(x)) > epsG & k < kmax
        %disp("H: " + h(x));
        [Q, La] = eig(h(x)); 
        % disp("Q: " + Q); % Elements de Q
        % disp("La: " + La); % Matriu Diagonal de Vaps
        La = diag(La); % Agafem només els elements de la diagonal
        if La(1:1) < 0 % Comprovem si hi ha VAPS negatius
            vaps_negatius = true;
        end

        La_prima = max(delta, La);
        % disp("La_prima: " + La_prima); % Vector de Vaps modificats
        La_Matrix = diag(La_prima);
        disp("La_Matrix: " + La_Matrix); % Matriu Diagonal de Vaps modificada
        
        B = Q * La_Matrix * Q';
        %disp("B: " + B); % Hessiana modificada
        d = - B^-1 * g(x);

        [al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW); 
       
        x = x + al*d; 
        Hk(:,:,end+1) = [B];
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al]; iWk = [iWk, iWout]; 
        k = k + 1;
    end
    if vaps_negatius == true
        disp(" --- !!!TENIM VAPS NEGATIUS, HESSIANA NO + DEF.!!! --- ");
    end
end