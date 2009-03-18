%
% Author:   Gary Mirams
% Date:     11/6/2005

close all
clear


% mitogenic_factorF = 5e-4;
mitogenic_factorF = 1.0/25.0;
mutation = 0;
wnt_zero = 1.0;
tstart = 0.0;
tfinal = 400;
hypothesis = 1;
n1=100;
n2=10000;
n3=1000;

%
% The variables are
%
% 1. r = pRb
% 2. e = E2F1 (This is the S-phase indicator)
% 3. i = CycD (inactive)
% 4. j = CycD (active)
% 5. p = pRb-p
% 6. D = APC destruction complex
% 7. X = Axin
% 8. Cu = Beta Cat marked for ubiquitination
% 9. Co = Open form Beta Cat
% 10. Cc = Closed form Beta Cat
% 11. Mo = Open form Mutant Beta Cat
% 12. Mc = Closed form Mutant Beta Cat
% 13. A = Free Adhesion molecules
% 14. Ca = BetaCat/Adhesion
% 15. Ma = Mutant BetaCat/Adhesion
% 16. T = free TCF
% 17. Cot = Open BetaCat/TCF
% 18. Cct = Closed BetaCat/TCF
% 19. Mot = Open Mutant BetaCat/TCF
% 20. Mct = Closed Mutant BetaCat/TCF
% 21. Y = Wnt Target protein
% 22. Wnt level
% Set up original parameters

k1 = 1.0;
k2 = 1.6;
k3 = 0.05;
k16 = 0.4;
k34 = 0.04;
k43 = 0.01;
k61 = 0.3;
k67 = 0.0; % were k67 = 0.7; k76 = 0.1;
k76 = 0.0; % Reduced to simplify equations
k23 = 0.3;
k25 = 0.9;
k28 = 0.06;
k89 = 0.07;
k98 = 0.01;
a = 0.04;
J11 = 0.5;
J12 = 5.0;
J15 = 0.001;
J18 = 0.6;
J61 = 5.0;
J62 = 8.0;
J65 = 6.0;
J68 = 7.0;
J13 = 0.002;
J63 = 2.0;
Km1 = 0.5;
Km2 = 4.0;
Km4 = 0.3;
Km9 = 0.005;
kp = 0.05;
phi_pRb = 0.005;
phi_E2F1 = 0.1;
phi_CycDi = 0.023;
phi_CycDa = 0.03;
phi_AP1 = 0.01;
phi_pRbp = 0.06;
phi_pRbpp = 0.04;
phi_CycEi = 0.06;
phi_CycEa = 0.05;

% Non-dimensionalise ....
k2d = k2/(Km2*phi_E2F1);
k3d = k3*mitogenic_factorF/(Km4*phi_E2F1);
k34d = k34/phi_E2F1;
k43d = k43/phi_E2F1;
k23d = k23*Km2/(Km4*phi_E2F1);
ad = a/Km2;
J11d = J11*phi_E2F1/k1;
J12d = J12*phi_E2F1/k1;
J13d = J13*phi_E2F1/k1;
J61d = J61*phi_E2F1/k1;
J62d = J62*phi_E2F1/k1;
J63d = J63*phi_E2F1/k1;
Km1d = Km1/Km2;
kpd = kp/(Km2*phi_E2F1);
phi_r = phi_pRb/phi_E2F1;
phi_i = phi_CycDi/phi_E2F1;
phi_j = phi_CycDa/phi_E2F1;
phi_p = phi_pRbp/phi_E2F1;
k16d = k16*Km4/phi_E2F1;
k61d = k61/phi_E2F1;

mSa = 20;   %  nM/h
mSca = 250; %  (nMh)^-1
mSc = 25;   %  nM/h
mSct = 30;  %  (nMh)^-1
mSd = 100;  %  h^-1
mSt = 10;   %  nM/h
mSx = 10;   %  nM/h
mSy = 10;   %  h^-1
mDa = 2;    %  h^-1
mDca = 350; %  h^-1
mDc = 1;    %  h^-1
mDct = 750; %  h^-1
mDd = 5;    %  h^-1
mDdx = 5;   %  h^-1
mDt = 0.4;  %  h^-1
mDu = 50;   %  h^-1
mDx = 100;  %  h^-1
mDy = 1;    %  h^-1
mKc = 200;  %  nM
mKd = 5;    %  nM
mKt = 50;   %  nM
mPc = 0.0;  %  h^-1
mPu = 100;  %  h^-1
mXiD = 5;   %  h^-1
mXiDx = 5;  %  h^-1
mXiX = 200; %  h^-1
if (hypothesis == 1)
    mXiC = 0.0; %  h^-1 (0.0 FOR HYPOTHESIS ONE, 5000 for Hypothesis two)
elseif (hypothesis ==2)
    mXiC = 5000.0;
else
    disp('Error');
end
parameters = [k2d k3d k34d k43d k23d ad J11d J12d J13d J61d J62d J63d ...
    Km1d kpd phi_r phi_i phi_j phi_p k16d k61d mutation mSa mSca mSc ...
    mSct mSd mSt mSx mSy mDa mDca mDc mDct mDd mDdx mDt mDu mDx mDy ...
    mKc mKd mKt mPc mPu mXiD mXiDx mXiX mXiC]';

wnt_zero_a = linspace(0,0.63,n1)';
wnt_zero_b = linspace(0.63000001,0.63999999,n2)';
wnt_zero_c = linspace(0.64,1,n3)';
wnt_zero = [wnt_zero_a; wnt_zero_b; wnt_zero_c];
n = n1+n2+n3;
for i=1:n
    d_d_hat = mDd + mXiD*wnt_zero(i);
    d_d_x_hat = mDdx + mXiDx*wnt_zero(i);
    d_x_hat = mDx + mXiX*wnt_zero(i);
    p_c_hat = mPc + mXiC*wnt_zero(i);

    sigma_D = 0.0;   % for healthy cells
    sigma_B = 0.0;   % for healthy cells

    % Mutations take effect by altering the level of beta-catenin
    if (abs(mutation - 0.0)<1e-6) % Healthy
    elseif (abs(mutation - 1.0)<1e-6) % APC +/-
        sigma_D = 0.5;
    elseif (abs(mutation - 2.0)<1e-6) % Beta-Cat D45
        sigma_B = 0.5;
    elseif (abs(mutation - 3.0)<1e-6) % APC -/-
        sigma_D = 1.0;
    else	% Shouldn't happen
        disp('Wnt Cell Cycle model using invalid mutation state');
    end

    steady_D = ((1.0-sigma_D)*mSd*mSx)/((1.0-sigma_D)*mSd*d_d_hat + ...
        d_x_hat*(d_d_hat + d_d_x_hat));
    steady_X = (mSx*(d_d_hat+d_d_x_hat))/((1.0-sigma_D)*mSd*d_d_hat+ ...
        d_x_hat*(d_d_hat+d_d_x_hat));
    steady_Cf = ((mSc-mDc*mKd - mPu*steady_D)+ ...
        sqrt((mSc-mDc*mKd - mPu*steady_D)^2 + (4.0*mSc*mDc*mKd)))/(2.0*mDc);
    steady_Cu = (mPu*steady_D*steady_Cf)/(mDu*(steady_Cf+mKd));
    theta = mDc+ (mPu*steady_D)/(steady_Cf + mKd);
    steady_Co = ( mSc - p_c_hat - theta*mKc + ...
        sqrt(4.0*mSc*theta*mKc + (mSc - p_c_hat - theta*mKc)^2) )/(2.0*theta);
    steady_Cc = steady_Cf - steady_Co;
    steady_Mo = 0.0;
    steady_Mc = 0.0;
    steady_A = mSa/mDa;
    steady_Ca = mSa*mSca*steady_Co/(mDa*mDca);
    steady_Ma = 0.0;
    steady_T = mSt/mDt;
    steady_Cot = mSct*mSt*steady_Co/(mDt*mDct);
    steady_Cct = mSct*mSt*steady_Cc/(mDt*mDct);
    steady_Mot = 0.0;
    steady_Mct = 0.0;
    steady_Y = (mSct*mSt*mSy*steady_Cf)/(mDy*(mSct*mSt*steady_Cf + mDct*mDt*mKt));

    % Set up Initial Conditions (from Figure 3 caption)
    pRb_zero = 7.357;
    E2F1_zero = 0.6852;
    CycDi_zero = 0.0207;
    CycDa_zero = 0.001;
    pRbp_zero = 0.001;

    r_zero = pRb_zero*phi_E2F1/k1;
    e_zero = E2F1_zero/Km2;
    i_zero = CycDi_zero/Km4;
    j_zero = CycDa_zero/Km4;
    p_zero = pRbp_zero*phi_E2F1/k1;

    xzero = [r_zero e_zero i_zero j_zero p_zero steady_D steady_X steady_Cu ...
        steady_Co steady_Cc steady_Mo steady_Mc steady_A steady_Ca steady_Ma ...
        steady_T steady_Cot steady_Cct steady_Mot steady_Mct steady_Y wnt_zero(i)]';


    options = odeset('Events',@stopping_condition,'AbsTol',1e-5,'RelTol',1e-3);
    % options = odeset('AbsTol',1e-5,'RelTol',1e-3);
    [t y] = ode15s(@IngeSwat2004ChasteEquations,[tstart tfinal],xzero,options,parameters);
    % Work out a transition time for this level of beta-catenin...
    % Detect when G1/S transition happens (when E2F1 > 1 ish)
    for j=1:length(t)
        if y(j,2) > 1.0
            transition_time = t(j-1)+(1.0-y(j-1,2))/...
                (y(j,2)-y(j-1,2))*(t(j)-t(j-1));
            break
        end
        % Catch if it doesn't reach transition
        if j==length(t)
            transition_time = 10000.0;
        end
    end
    transition_result(i) = transition_time + 8.5; % time for other phases.
end

crypt_height = (30.0*sqrt(3)/2.0)*(20.1/23.0);
crypt_height_cutoff = (8.0/18.0)*crypt_height;

top_of_crypt = 23;

figure
subplot(3,1,1)
height_zero = (1.0 - wnt_zero).*(crypt_height_cutoff);
plot(height_zero,transition_result,'k.')
axis([0 top_of_crypt 14 26])
set(gca,'Fontsize',14)
xlabel('Height (cells)','Fontsize',16)
ylabel('Cell cycle duration (hours)','Fontsize',16)

% Plot simulation results
data = load('wnt_vs_cycle.dat');
subplot(3,1,2)
scatter(data(:,1),data(:,3),'kx')
set(gca,'Fontsize',14)
axis([0 top_of_crypt 14 26])
xlabel('Height (cells)','Fontsize',16)
ylabel('Cell cycle duration (hours)','Fontsize',16)

% PLot the Wnt Concentration

n=1000;
height = linspace(0,top_of_crypt,n);
for i=1:n
    if height(i) < crypt_height_cutoff
        wnt(i) = 1.0 - height(i)/crypt_height_cutoff;
    else
        wnt(i) = 0.0;
    end
end
subplot(3,1,3)
plot(height,wnt,'k-','LineWidth',2)
axis([0 top_of_crypt 0 1])
set(gca,'Fontsize',14)
xlabel('Height (cells)','Fontsize',16)
ylabel('Wnt Concentration','Fontsize',16)

% figure
% subplot(2,1,1)
% plot(t,y(:,2),'r-', t,y(:,1),'b--')
% xlabel('time (hours)')
% ylabel('E2F1 (red)')
% subplot(2,1,2)
% plot(t,y(:,17)+y(:,18),'r-', t,y(:,19)+y(:,20),'b--')
% xlabel('time (hours)')
% ylabel('beta-cat (red), beta-cat mutant (blue)')
