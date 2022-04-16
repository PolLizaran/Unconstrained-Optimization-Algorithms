function [B, tau] = B_MODIF_CMI(x, h)
    tau = 0; taus =[tau]; it = 0; n = size(x); n1 = n(1);
    La_UB = norm(h(x),'fro');
    disp("La_UB " + La_UB) % Cota superior de lambda_n
    exists_RtR = inf;
    while exists_RtR > 0
        [R, exists_RtR] = chol(h(x) + tau*eye(n1));
        it = it + 1;
        tau = (1.01 - 1/(2^it))*La_UB; taus = [taus, tau];
    end
    tau = taus(end - 1);
    disp("tau " + tau) % Tau més petita per la qual s'ha pogut factoritzar
    disp("R " + R) % Factor de Cholesky
    B = R'*R; 
    disp("B " + B) % Hessiana modificada
end
