%
% Inge's ODEs coupled to the Swat equations.
%
function dx=IngeSwat2004ChasteEquations(time,x,parameters)
% 
% The variables are
% 1. pRb
% 2. E2F1
% 3. CycD (inactive)
% 4. CycD (active)
% 5. pRb-p
% 6. D
% 7. X
% 8. Cu
% 9. Co
% 10. Cc
% 11. Mo
% 12. Mc
% 13. A
% 14. Ca
% 15. Ma
% 16. T
% 17. Cot
% 18. Cct
% 19. Mot
% 20. Mct
% 21. Y
% 22. Wnt level
%
dx = zeros(22,1);

k2 = parameters(1);
k3 = parameters(2);
k34 = parameters(3);
k43 = parameters(4);
k23 = parameters(5);
a = parameters(6);
J11 = parameters(7);
J12 = parameters(8);
J13 = parameters(9);
J61 = parameters(10);
J62 = parameters(11);
J63 = parameters(12);
Km1 = parameters(13);
kp = parameters(14);
phi_r = parameters(15);
phi_i = parameters(16);
phi_j = parameters(17);
phi_p = parameters(18);
k16 = parameters(19);
k61 = parameters(20);

% Deal with different mutations
mutation = parameters(21);

% Initialize Inge's model parameters
mSa = parameters(22);
mSca = parameters(23);
mSc = parameters(24);
mSct = parameters(25);
mSd = parameters(26);
mSt = parameters(27);
mSx = parameters(28);
mSy = parameters(29);
mDa = parameters(30);
mDca = parameters(31);
mDc = parameters(32);
mDct = parameters(33);
mDd = parameters(34);
mDdx = parameters(35);
mDt = parameters(36);
mDu = parameters(37);
mDx = parameters(38);
mDy = parameters(39);
mKc = parameters(40);
mKd = parameters(41);
mKt = parameters(42);
mPc = parameters(43);
mPu =parameters(44); 
mXiD =parameters(45);
mXiDx = parameters(46);
mXiX = parameters(47);
mXiC = parameters(48);

r = x(1);
e = x(2);
i = x(3);
j = x(4);
p = x(5);
D = x(6);
X = x(7);
Cu = x(8);
Co = x(9);
Cc = x(10);
Mo = x(11);
Mc = x(12);
A = x(13);
Ca = x(14);
Ma = x(15);
T = x(16);
Cot = x(17);
Cct = x(18);
Mot = x(19);
Mct = x(20);
Y = x(21);
WntLevel= x(22);

% mXiC = 5000; % This line makes the equations use hypothesis two
% 
% if (time <= 1)
%     WntLevel = 1.0 - time;
% else
%     WntLevel = 0.0;
% end

Cf = Cc+Co;
Ct = Cct+Cot;
Mf = Mc+Mo;
Mt = Mct+Mot;
    
d_d_hat = mDd + mXiD*WntLevel;
d_d_x_hat = mDdx + mXiDx*WntLevel;
d_x_hat = mDx + mXiX*WntLevel;
p_c_hat = mPc + mXiC*WntLevel;
    
sigma_D = 0.0;   % for healthy cells
sigma_B = 0.0;   % for healthy cells

% Mutations take effect by altering the level of beta-catenin
if (abs(mutation - 0.0)<1e-6) % healthy
    % DO Nowt.
elseif (abs(mutation - 1.0)<1e-6) % APC +/-
    sigma_D = 0.5;
elseif (abs(mutation - 2.0)<1e-6) % Beta-Cat D45
    sigma_B = 0.5;
elseif (abs(mutation - 3.0)<1e-6) % APC -/-
    sigma_D = 1.0;
else	% Shouldn't happen
    disp('Wnt Cell Cycle model using invalid mutation state');
end

% dr
dx(1) = e/(Km1+e)*J11/(J11+r)*J61/(J61+p) - k16*r*j+k61*p-phi_r*r;
% de
dx(2) = kp+k2*(a^2+e^2)/(1+e^2)*J12/(J12+r)*J62/(J62+p)-e;
% di
dx(3) = k3*(Ct+Mt)+k23*e*J13/(J13+r)*J63/(J63+p)+k43*j-k34*i*j/(1+j)-phi_i*i;
% dj
dx(4) = k34*i*j/(1+j) - (k43+phi_j)*j;
% dp
dx(5) = k16*r*j-k61*p-phi_p*p;

% Inge's Wnt ODEs
dx(6) = (1.0-sigma_D)*mSd*X - (d_d_hat + d_d_x_hat)*D;
dx(7) = mSx - (1.0-sigma_D)*mSd*X - d_x_hat*X + d_d_x_hat*D;
dx(8) = (mPu*D*Cf)/(Cf+mKd) - mDu*Cu;
dx(9) = (1.0-sigma_B)*mSc + mDca*Ca + mDct*Cot - (mSca*A + mSct*T + mDc)*Co ...
        - (p_c_hat*Co)/(Co + Mo + mKc) - (mPu*D*Co)/(Cf+mKd);
dx(10) = (p_c_hat*Co)/(Co + Mo + mKc) + mDct*Cct - (mSct*T + mDc)*Cc ...
        - (mPu*D*Cc)/(Cf+mKd);
dx(11) = sigma_B*mSc + mDca*Ma + mDct*Mot - (mSca*A + mSct*T + mDc)*Mo...
        - (p_c_hat*Mo)/(Co + Mo + mKc);
dx(12) = (p_c_hat*Mo)/(Co + Mo + mKc) + mDct*Mct - (mSct*T + mDc)*Mc;
dx(13) = mSa + mDca*(Ca+Ma) - (mSca*(Co+Mo) + mDa)*A;
dx(14) = mSca*Co*A - mDca*Ca;
dx(15) = mSca*Mo*A - mDca*Ma;
dx(16) = mSt + mDct*(Ct+Mt) - mSct*(Cf+Mf)*T - mDt*T;
dx(17) = mSct*Co*T - mDct*Cot;
dx(18) = mSct*Cc*T - mDct*Cct;
dx(19) = mSct*Mo*T - mDct*Mot;
dx(20) = mSct*Mc*T - mDct*Mct;
dx(21) = (mSy*(Ct+Mt))/(Ct + Mt + mKt) - mDy*Y;

dx(22) = 0.0;   % don't change Wnt level

dx(1) = dx(1)*(0.1*60.0);   % make into d/dt in hours.
dx(2) = dx(2)*(0.1*60.0);
dx(3) = dx(3)*(0.1*60.0);
dx(4) = dx(4)*(0.1*60.0);
dx(5) = dx(5)*(0.1*60.0);
