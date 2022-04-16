% Funció que retorna la longitud de pas que ha de seguir l'algoritme i el
% paràmetre iWout que determina quina de les condicions de Wolfe compleix

function [al, iWout] = BLS(x, d, f, g, h, almax, almin, rho, c1, c2, iW)
    if iW == 0 % ELS 
        al = -(g(x)'*d)/(d'*h(x)*d); 
        iWout = 5; % No comprovem les condicions de Wolfe
    else 
        al = almax; % Inicialització de la alpha
        [satisfy, iWout] = WOLFE(x, al, d, f, g, c1, c2, iW);
        while ~satisfy & al >= almin % Hem de continuar
            al = rho*al; % Decrementem alpha amb un factor 'rho'
            [satisfy, iWout] = WOLFE(x, al, d, f, g, c1, c2, iW);
        end
    end
end