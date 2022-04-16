%﻿ iW= 0: exact LS; IW=1: BLS with WC; iW=2: BLS with SWC
% Hk = Bk per isd = 5,6

function [xk,dk,alk,iWk,betak,Hk,tauk]= ...
          uo_solve(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW, isd, icg, irc, nu, delta)
    fprintf(' \n Inici  [uo_solve]\n'); 
    al_meth = ["ELS", "BLS", "BLS", "BLS"];
    disp("Alpha es calcula amb " + al_meth(iW + 1));
    xk = []; dk = []; alk = []; iWk = []; betak = []; Hk = []; tauk = [];
    if isd == 1 % GM 
        fprintf('Apliquem GM \n');
        [xk, dk, alk, iWk] = ...
          GM(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW);
    elseif isd == 2 % CGM
        fprintf('Apliquem CGM amb condició de reinici RC%d, i càlcul de la beta amb %d \n', irc, icg);
        [xk, dk, alk, iWk, betak] = ...
          CGM(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW, icg, irc, nu);
    elseif isd == 3 % BFGS
        fprintf('Apliquem BFGS \n');
        [xk, dk, alk, iWk, betak, Hk] = ...
            BFGS(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW);
    elseif isd == 4 % NM
        fprintf('Apliquem NM \n');
        [xk, dk, alk, iWk, Hk] = ...
            NM(x, g, h, epsG, kmax);
    elseif isd == 5 % MNM-SD 
        fprintf('Apliquem MNM-SD amd delta = %d \n', delta);
        [xk, dk, alk, iWk, Hk] = ...
            MNM_SD(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW, delta); 
    elseif isd == 6 % MNM-CMI
        fprintf('Apliquem MNM-CMI \n');
        [xk, dk, alk, iWk, Hk, tauk] = ...
            MNM_CMI(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW); 
    end
    fprintf(' \n Final  [uo_solve]\n\n');
end



%[xk, dk, alk, iWk] = GM_BLS(x, f, g, h, epsG, kmax, almax, almin, rho, c1, c2, iW);