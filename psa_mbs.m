function price = psa_method()
    r0 = 0.078; r_bar = 0.08; k = 0.6; sig = 0.12; T = 30;
    PV = 100000; WAC = 0.08;
    psa_method(r0, r_bar,k  , sig, T, PV, WAC)
end

function price = psa_method(r0, r_bar, k , sig, T, PV, WAC)
    rng(1)
    RH = 40; % Rate horizon (in years)
    dt = 1/12;
    N = RH/dt; M = 1000; n = T*12; r = WAC/12;
    z = normrnd(0,1,M, N)*sqrt(dt);
    rt = [];
    rt(1:M,1) = r0;
    for(i = 1:N)
        rt(:, i+1) = rt(:, i) + k.*(r_bar - rt(:,i) ).*dt + sig.*sqrt(rt(:,i)).*z(:,i);
    end
    prices =[];
    for(i = 1:M)
        cfs = [];
        PV_t_1 = PV;
        for(t = 1:30*12)
            PMT = PV_t_1 * r / ( 1- 1/(1+r)^n );
            CPR_t = cpr_psa(t);
            ct = PV_t_1*r/ (1-(1+r)^(-n+(t-1))) ...
               + (PV_t_1 - PV_t_1*r*(1/(1-(1+r)^(-n + t-1)) - 1 ) )*(1-(1-CPR_t)^(1/12));
            tpp = PV_t_1*r * (1 / (1-(1+r)^(-n+(t-1))) - 1) ...
               + (PV_t_1 - PV_t_1*r*(1/(1-(1+r)^(-n + t-1)) - 1 ) )*(1-(1-CPR_t)^(1/12));
            cfs(t) = exp(-sum(rt(i,1:t)) * dt) *ct;
            PV_t_1  = PV_t_1 - tpp;
        end
        prices(end+1) = sum(cfs);
    end
    price = mean(prices);
end


function cpr = cpr_psa(t)
    if(t<=30)
        cpr=0.002*t;
    else
        cpr = 0.06;
    end
end