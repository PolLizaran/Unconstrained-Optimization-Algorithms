% Donat un punt retorna el coeficient de beta tenint en compte condicions
% d'aturada i el dos mètodes estudiats per la tria de beta.

function [beta] = BETA_CHOICE(x, xk, g, icg, irc, nu, k)
    % Stopping conditions (irc=0 -> (no restart); irc=1 -> RC1; irc=2 -> RC2)
    restarting = false;
    if irc == 1 & mod(k, length(x)) == 0
        restarting = true; beta = 0;
    elseif irc == 2 & (abs(g(x)'*g(xk(:,end)))/norm(g(x))^2 >= nu)
        restarting = true; beta = 0;
    end
    
    % Choice of the beta coeficients
    if ~restarting & icg == 1 % FR (Fletcher-Reeves)
        beta =  (g(x)'*g(x))/norm(g(xk(:,end)))^2;
    elseif ~restarting & icg == 2 % PR (Polak-Ribière)
        beta = max(0 , g(x)'* (g(x) - g(xk(:,end)))/( norm(g(xk(:,end)))^2 ));
    end
end

