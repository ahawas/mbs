function price = numerix_mbs()
    r0 = 0.078; r_bar = 0.08; k = 0.6; sig = 0.12; T = 30;
    PV = 100000; WAC = 0.08;
    numerix_method(r0, r_bar,k  , sig, T, PV, WAC)
end

function price = numerix_method(r0, r_bar, k , sig, T, PV, WAC)
    RH = 40; % Rate horizon (in years)
    dt = 1/12;
    N = RH/dt; M = 1000; n = T*12; r = WAC/12;
    z = normrnd(0,1,M/2, N)*sqrt(dt);
    z= [z ;-z];
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
            CPR_t = cpr_numerix(t, rt, dt, i, PV_t_1, PMT, r, sig, k, r_bar);
            %CPR_t = cpr_psa(t);      
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

function cpr=cpr_numerix(t, rt, dt, int_path, PV_t_1, PMT,r, sig, k, r_bar)
    n = 30*12;
    PV_0 = 100000;
    BU_t = 0.3 + 0.7* PV_t_1 / PV_0;
    SG_t = min(1,t/30);
    R = r*12;
    months = [0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.1, 1.18, 1.22, 1.23, 0.98];
    if(mod(t,12) == 0)
        SY_t = months(12);
    else
        SY_t = months(mod(t,12));
    end
    %r_t_1_10 = 1/10 * sum(rt(int_path,t:t+10*12)*dt);
    r_t_1_10 = 1/10 * -log(cir_exp(0,10, rt(int_path,t), r_bar, sig, k));
    RI_t = 0.28 + 0.14*atan(-8.57 + 430*(R - r_t_1_10));
    cpr = RI_t * BU_t*SG_t*SY_t;
end
function price = cir_exp(t,T, rt, r_bar, sig, k )
    theta = (k^2 + 2*sig^2)^0.5;
    phi = 2*theta / (sig^2*(exp(theta*(T-t))-1)) ;
    psi = (k + theta) / sig^2;
    h1 = theta;h2 = (k+h1)/2; h3 = 2*k*r_bar/(sig^2);
    A_T_S = (h1*exp(h2*(T-t) )/...
         (h2*(exp(h1*(T-t)) -1 ) + h1) )^h3;
    B_T_S= (exp(h1*(T-t))-1)/...
         (h2*(exp(h1*(T-t))-1)+h1) ;
    price = A_T_S * exp(-B_T_S*rt);
end