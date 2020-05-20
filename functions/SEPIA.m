function [x,V]=SEPIA(PAR,V)

%rename parameters
n=V.n;  %number of nodes
deltaE  = 1/PAR(2);
deltaP  = 1/PAR(3);
sigma   = PAR(4);
eta     = 1/PAR(5);
gammaI  = 1/PAR(6);
alphaI  = 1/PAR(7);
alphaH  = 1/PAR(7);
gammaH  = 1/PAR(6);
epsilonA= PAR(8);  %for brevity the ratio \beta_A/\beta_P is termed epsilonA
r       = PAR(9);  %the model assumes that clesses S E P A R have the same r values
Deltat0 = PAR(10);
epsilonI= PAR(11)*epsilonA; %for brevity the ratio \beta_I/\beta_P is termed epsilonI
betaP1P0= PAR(13);
betaP2P1= PAR(14);
gammaQ  = V.gammaQ_over_gammaH*gammaH;
gammaA  = V.gammaA_over_gammaQ*gammaQ;

% Parameter betaP0 is expressed as a function of the local reproductive
% number R_0 and the other relevnt parameters
betaP0 = PAR(1)/(1/deltaP + epsilonI*sigma/(gammaI + alphaI + eta) + epsilonA*(1-sigma)/gammaA);

% Add initial conditions in the exposed class for the seeding nodes
V.x0(n+V.seeding) = 10.^(PAR(V.nPAR_model+V.n_reg+1:V.nPAR_model+V.n_reg+length(V.seeding)));

% Define simulation time (time_model)
V.time_model = (V.Date(1)-(Deltat0)):V.time_model_final;

% Calculate mobility ratio for each node (1st dimension) and day of
% V.time_model (second dimension) through linear interpolation of V.mob
mob_ratio = interp1(V.tmob, V.mob',V.time_model)';

% Calculate  trasmission ratio (beta_ratio) for each node (1st dimension) and day of
% V.time_model (second dimension) through linear interpolation of V.mob
on = ones(V.n,1);
betaP3P2_reg = PAR(V.nPAR_model+1:V.nPAR_model+V.n_reg)'; %regional value
betaP3P2 = betaP3P2_reg(V.prov_IDreg); %province value
beta_p = eval(V.beta_string);
beta_ratio = interp1(V.tbeta,beta_p',V.time_model)';

% run simulation
[t,x] = ode45(@eqs,V.time_model,V.x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function with model equations
    function dxdt=eqs(t,x)
        i_t = floor(t-V.time_model(1))+1;  %index of time (day of simulation)
        
        %compartments
        S = x(1:n);       %susceptible
        E = x(n+1:2*n);   %exposed
        P = x(2*n+1:3 *n);%peak infectivity
        I = x(3*n+1:4*n); %heavily symptomatic
        A = x(4*n+1:5*n); %mildly symptomatic/asymptomatic
        Q = x(5*n+1:6*n); %quarantined
        H = x(6*n+1:7*n); %hospitalized
        R = x(7*n+1:8*n); %recovered
        
        % Compute mobility matrix C (see model equations). According to the
        % assumptions used in the preprint: C_S=C_E=C_P=C_A=C_R. Moreover,
        % C_I is an idenity matrix (no extra node mobility for the I class)
        C = (r*V.p.*mob_ratio(:,i_t)).*V.q; %compute off diagonal elements
        C(1:n+1:end) = 1-(sum(C,2)-diag(C));%set diagonal imposing row-stocasticity
        
        % Force of infection lambda (this formulation is a numerically faster simplification of
        % the equation reported in the preprint text valid under the
        % assumptions listed above for the computation of C
        lambda=C*((C'*(betaP0*beta_ratio(:,i_t).*(P+epsilonA*A))+epsilonI*betaP0*beta_ratio(:,i_t).*I)./(C'*(S+E+P+R+A)+I));
        
        % Model equations
        dSdt = -lambda.*S;
        dEdt = lambda.*S-deltaE*E;
        dPdt = deltaE*E-deltaP*P;
        dIdt = sigma*deltaP*P-(eta+gammaI+alphaI)*I;
        dAdt = (1-sigma)*deltaP*P-gammaA*A;
        dQdt = V.zeta*eta*I-gammaQ.*Q;
        dHdt = (1-V.zeta)*eta*I-(gammaH+alphaH)*H;
        dRdt = gammaI*I+gammaA*A+gammaH*H+gammaQ*Q;
        dRRdt= gammaH*H;                            % recovered from the hospital
        dDdt = alphaH*H;                            % total death
        dcumHdt = (1-V.zeta)*eta*I;                 % cumulative hospitalized cases
        
        dxdt = [dSdt; dEdt; dPdt; dIdt; dAdt; dQdt; dHdt; dRdt; dRRdt; dDdt; dcumHdt];
        return
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end