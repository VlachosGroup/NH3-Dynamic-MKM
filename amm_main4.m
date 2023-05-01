function [tr, sr, Conv] = amm_main4(T_in)
%%          -------------------------------------------------
%                        NH3  Micro-kinetic model
%                         Vlachos Research Group
%                 Chemical and Biomolecular Egineering
%                         University of Delaware
%
%             Gerhard R Wittreich, P.E.  (February 10, 2017)
%           --------------------------------------------------
%
%  Main program:
%      requires: ammonia.m      : Reaction ODE's and dynamical equations
%                amm_kinetics.m : Provides forward and reverse rate constants
%                amm_thermo.m   : NASA polynomials provides enthalpy and entropy
%                amm_coverage.m : Calculates lateral interactions of surface species
%                amm_BEP_LSR.m  : Adjusts thermo data for target metal
%                                 catalyst, provides Ea and omega for each reaction
%%
clearvars -except T_in
close all
fprintf ('Temperature = %3d\n',T_in)
%  Set key model parameters
%
global T T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 abyv...
    c_tot Stoic_surf Stoic_gas Stoic MWON Isobaric Ea A Stick R_e ...
    R_k R MW_N2 MW_H2 MW_NH3 SDEN_T SDEN_S SDTOT RR Q_name period...
    strain strain_pulse tspan Strain_Coef_H Strain_Coef_S A_LSR t_check tspan_cum...
	RR_Count RR_All tt
T_orig = T_in;                   % Reactor bulk temperature [K]
T = T_orig;                     % Initial catalyst temperature [K]
T_gas = T_orig;                 % Initial gas temperature [K]
P = 50.0;                        % Reactor pressure [atm]
beta = [0 1 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 1 1 1 1 1]';
Ea   = zeros(22,1);
Stick= [0.5 0.5 0.5 0.5 0.5 0.5]';
period = 1./2000;               % Pulse period [sec]
MWON = 0;                       % Motz-Wise Correction: 0 = Off (S); 1 = On (S/(1-S/2))
Isobaric = 1;                   % Isobaric Reactor (0=Isochoric, 1=Isobaric)
t_check = 0;
MW_H = 1.00797;                 % ---
MW_N = 14.0067;                 %
MW_H2 = 2*MW_H;                 % Molecular mass [g/mol]
MW_N2 = 2*MW_N;                 %
MW_NH3 = MW_N + 3*MW_H;         % ---
X_H2  = 3;                      %
X_N2  = 1;                      % Initial mole fraction in reactor feed
X_NH3 = 0;                      %
Y_H2  = X_H2 /(X_H2+X_N2+X_NH3);% Mole fractions
Y_N2  = X_N2 /(X_H2+X_N2+X_NH3);% normalized
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);% to 1
abyv = 1200.;                    % Catalyst loading (cm2 catalyst/cm3 reac volume)
V = 1.0;                        % Reactor volume (cm3)
%Q_in = 1.0*T_orig/298.15/P;    % 0 = Batch Reactor,  Any other value = CSTR [cm3/s]
Q_in = 1.0;                    % 0 = Batch Reactor,  Any other value = CSTR [cm3/s]
% SDEN_T = 2.16715e-09;           % Catalyst terrace site density (moles/cm2)
% SDEN_S = 4.4385e-10;            % Catalyst step site density (moles/cm2)
% SDTOT = SDEN_T + SDEN_S;        % Total catalyst site density (moles/cm2)

SDTOT = 2.6188e-9;              % Total catalyst site density (moles/cm2)
RATIO_S = 0.10;                 % Fraction of Step sites
SDEN_T = (1-RATIO_S) * SDTOT;   % Catalyst terrace site density
SDEN_S = RATIO_S * SDTOT;       % Catalyst step site density


i_strain = 0.00;                % Initial catalyst lattice strain
R_e = 1.987e-3;                 % Gas constant, (kcal/mol K)
R_k = 8.31451e7;                % Gas constant, (g cm2/mol K s)
R = 82.057;                     % Gas constant, (cm3 atm/K mol)
c_tot = P/(R*T);                % Total starting moles in reactor [mol/cm3]
c_H2 = Y_H2*c_tot;              % Moles H2 in reactor [mol/cm3]
c_N2 = Y_N2*c_tot;              % Moles N2 in reactor [mol/cm3]
c_NH3 = Y_NH3*c_tot;            % Moles NH3 in reactor [mol/cm3]
A    = [0 0 0 0 0 0 0 0 1.56E19 1.56E19 1.56E19 1.56E19 1.56E19 0 0 1.56E19]';
    
% Strain impact on enthalpy [kcal/mol]
Strain_H_N2_T =  [ 3.35238230  0.  -4.74700358]/R_e;
Strain_H_N_T =   [ 12.3535236  0. -11.07668235]/R_e;
Strain_H_H_T =   [ 2.65958231  0.  -2.39723156]/R_e;
Strain_H_NH3_T = [ 2.05050749  0.  -4.03832726]/R_e;
Strain_H_NH2_T = [ 2.90623934  0.  -3.51266083]/R_e;
Strain_H_NH_T =  [ 2.69136079  0.  -6.38940051]/R_e;
Strain_H_N2_S =  [ 4.74665949  0.  -2.74674134]/R_e;
Strain_H_N_S =   [ 2.55994332  0.  -2.99094398]/R_e;
Strain_H_N_S3 =  [ 6.07410400  0.  -6.97056453]/R_e;
Strain_H_H_S =   [-0.58804601  0.  -0.00586277]/R_e;
Strain_H_NH3_S = [ 0.42259649  0.  -0.35136283]/R_e;
Strain_H_NH2_S = [-1.43738869  0.   0.66868959]/R_e;
Strain_H_NH_S =  [ 0.01759582  0.  -1.18782454]/R_e;

% Strain impact on entropy [kcal/mol/K]
Strain_S_N2_T =  [-1.25416455  0.  -0.79706019]/R_e/1000;
Strain_S_N_T =   [ 0.22363829  0.   0.04043008]/R_e/1000;
Strain_S_H_T =   [-0.11958374  0.  -0.03086301]/R_e/1000;
Strain_S_NH3_T = [-0.77580104  0.  -0.75740069]/R_e/1000;
Strain_S_NH2_T = [-0.55928909  0.   0.11301182]/R_e/1000;
Strain_S_NH_T =  [-0.38660173  0.   0.09730831]/R_e/1000;
Strain_S_N2_S =  [ 0.19807204  0.   0.20206281]/R_e/1000;
Strain_S_N_S =   [ 0.00916901  0.   0.13214497]/R_e/1000;
Strain_S_N_S3 =  [ 0.89917975  0.  -0.13420608]/R_e/1000;
Strain_S_H_S =   [-0.34809169  0.   0.19211268]/R_e/1000;
Strain_S_NH3_S = [ 3.56664371  0.   1.40560633]/R_e/1000;
Strain_S_NH2_S = [-0.35706660  0.  -0.11538387]/R_e/1000;
Strain_S_NH_S =  [-0.27169163  0.   0.02404511]/R_e/1000;

Stoic_surf = [ 1  0  0  0  0  0  0  0  0 -1  0;... % Reaction
              -1  2  0  0  0  0  0  0  0 -1  0;... %
               0  0  2  0  0  0  0  0  0 -2  0;... % Surface
               0  0  1 -1  1  0  0  0  0 -1  0;... %
               0  0  1  0 -1  1  0  0  0 -1  0;... % Stoichiometry
               0  1  1  0  0 -1  0  0  0 -1  0;... %
               0  0  0  1  0  0  0  0  0 -1  0];   % ---
Stoic_surf = [Stoic_surf zeros(7,7)];
Stoic_surf = [Stoic_surf; zeros(7,10) Stoic_surf(:,1:6) Stoic_surf(:,10) zeros(7,1)];
Stoic_surf = [Stoic_surf;...
               0 -1  0  0  0  0  0  0  0  1  0  1  0  0  0  0 -1  0;...
               0  0 -1  0  0  0  0  0  0  1  0  0  1  0  0  0 -1  0;...
               0  0  0 -1  0  0  0  0  0  1  0  0  0  1  0  0 -1  0;...
               0  0  0  0 -1  0  0  0  0  1  0  0  0  0  1  0 -1  0;...
               0  0  0  0  0 -1  0  0  0  1  0  0  0  0  0  1 -1  0;...
               0  1  0  0  0  0  0  0  0 -1 -1  1  0  0  0  0  0  0;...
               0  0  0  0  0  0  0  0  0 -1 -1  1  0  0  0  0  0  1;...
               0 -1  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1  1];
Stoic_gas =  [ 0  0  0  0  0  0 -1  0  0  0;... % Reaction
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0 -1  0  0;... % Gas
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0  0  0  0;... % Stoichiometry
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0  0 -1  0];   % ---
Stoic_gas =  [Stoic_gas zeros(7,8); Stoic_gas zeros(7,8); zeros(8,18)];
Stoic = Stoic_surf + Stoic_gas;                 % Total stoichiometry
StrainCoef = [-0.04 0.0 0.04];
Strain_Coef_H = zeros(13,2);
    Strain_Coef_H(1,:)  = polyfit(StrainCoef,Strain_H_N2_T,1);
    Strain_Coef_H(2,:)  = polyfit(StrainCoef,Strain_H_N_T,1);
    Strain_Coef_H(3,:)  = polyfit(StrainCoef,Strain_H_H_T,1);
    Strain_Coef_H(4,:)  = polyfit(StrainCoef,Strain_H_NH3_T,1);
    Strain_Coef_H(5,:)  = polyfit(StrainCoef,Strain_H_NH2_T,1);
    Strain_Coef_H(6,:)  = polyfit(StrainCoef,Strain_H_NH_T,1);
    Strain_Coef_H(7,:)  = polyfit(StrainCoef,Strain_H_N2_S,1);
    Strain_Coef_H(8,:)  = polyfit(StrainCoef,Strain_H_N_S,1);
    Strain_Coef_H(9,:)  = polyfit(StrainCoef,Strain_H_H_S,1);
    Strain_Coef_H(10,:) = polyfit(StrainCoef,Strain_H_NH3_S,1);
    Strain_Coef_H(11,:) = polyfit(StrainCoef,Strain_H_NH2_S,1);
    Strain_Coef_H(12,:) = polyfit(StrainCoef,Strain_H_NH_S,1);
    Strain_Coef_H(13,:) = polyfit(StrainCoef,Strain_H_N_S3,1);

Strain_Coef_S = zeros(13,3);
    Strain_Coef_S(1,:)  = polyfit(StrainCoef,Strain_S_N2_T,2);
    Strain_Coef_S(2,:)  = polyfit(StrainCoef,Strain_S_N_T,2);
    Strain_Coef_S(3,:)  = polyfit(StrainCoef,Strain_S_H_T,2);
    Strain_Coef_S(4,:)  = polyfit(StrainCoef,Strain_S_NH3_T,2);
    Strain_Coef_S(5,:)  = polyfit(StrainCoef,Strain_S_NH2_T,2);
    Strain_Coef_S(6,:)  = polyfit(StrainCoef,Strain_S_NH_T,2);
    Strain_Coef_S(7,:)  = polyfit(StrainCoef,Strain_S_N2_S,2);
    Strain_Coef_S(8,:)  = polyfit(StrainCoef,Strain_S_N_S,2);
    Strain_Coef_S(9,:)  = polyfit(StrainCoef,Strain_S_H_S,2);
    Strain_Coef_S(10,:) = polyfit(StrainCoef,Strain_S_NH3_S,2);
    Strain_Coef_S(11,:) = polyfit(StrainCoef,Strain_S_NH2_S,2);
    Strain_Coef_S(12,:) = polyfit(StrainCoef,Strain_S_NH_S,2);
    Strain_Coef_S(13,:) = polyfit(StrainCoef,Strain_S_N_S3,2);
    
    A_Strain_N2_T  = [9.05e+17, 4.32e+17, 5.94e+17]*(1-0.17)/(1-RATIO_S);
    A_Strain_NH3_T = [1.32e+18, 1.08e+18, 1.26e+18]*(1-0.17)/(1-RATIO_S);
    A_Strain_NH2_T = [7.68e+18, 5.25e+18, 4.74e+18]*(1-0.17)/(1-RATIO_S);
    A_Strain_NH_T  = [3.58e+19, 1.53e+19, 1.17e+19]*(1-0.17)/(1-RATIO_S);
    A_Strain_N2_S  = [2.40e+19, 2.59e+19, 2.65e+19]*0.17/RATIO_S;
    A_Strain_NH3_S = [1.17e+18, 7.46e+18, 3.21e+18]*0.17/RATIO_S;
    A_Strain_NH2_S = [4.49e+19, 5.30e+19, 5.43e+19]*0.17/RATIO_S;
    A_Strain_NH_S  = [4.70e+19, 4.28e+19, 4.18e+19]*0.17/RATIO_S;
    A_Strain_N2_ST = [8.18e+18, 8.81e+18, 9.02e+18]*((1-0.17)+(0.17))/((1-RATIO_S)+RATIO_S);
    A_Strain_N2_S3 = [8.18e+18, 8.81e+18, 9.02e+18]*((1-0.17)+(0.17))/((1-RATIO_S)+RATIO_S);

A_LSR = zeros(10,3);
    A_LSR(1,:) = polyfit(StrainCoef,A_Strain_N2_T,2);  % N2(T) + RU(T) = TS4_N2(T) = 2N(T) + RU(B)
    A_LSR(2,:) = polyfit(StrainCoef,A_Strain_NH3_T,2); % NH3(T) + RU(T) = TS1_NH3(T) = H(T) + NH2(T) + RU(B)
    A_LSR(3,:) = polyfit(StrainCoef,A_Strain_NH2_T,2); % NH2(T) + RU(T) = TS2_NH2(T) = H(T) + NH(T) + RU(B)
    A_LSR(4,:) = polyfit(StrainCoef,A_Strain_NH_T,2);  % NH(T) + RU(T) = TS3_NH(T) = N(T) + H(T) + RU(B)
    A_LSR(5,:) = polyfit(StrainCoef,A_Strain_N2_S,2);  % N2(S) + RU(S) = TS4_N2(S) = 2N(S) + RU(B)
    A_LSR(6,:) = polyfit(StrainCoef,A_Strain_NH3_S,2); % NH3(S) + RU(S) = TS1_NH3(S) = H(S) + NH2(S) + RU(B)
    A_LSR(7,:) = polyfit(StrainCoef,A_Strain_NH2_S,2); % NH2(S) + RU(S) = TS2_NH2(S) = H(S) + NH(S) + RU(B)
    A_LSR(8,:) = polyfit(StrainCoef,A_Strain_NH_S,2);  % NH(S) + RU(S) = TS3_NH(S) = N(S) + H(S) + RU(B)
    A_LSR(9,:) = polyfit(StrainCoef,A_Strain_N2_ST,2); % N2(S) + RU(T) = TS4_N2(S) = N(S) + N(T) + RU(B)
    A_LSR(10,:) = polyfit(StrainCoef,A_Strain_N2_S3,2);% N2(S) + RU(T) = TS4_N2(S) = N(S) + N(S3) + RU(B)

% ODE Solver options
MaxStep0 = period/10;
options0 = odeset ('MaxStep',MaxStep0,'NonNegative',...
    [1 2 3 4 5 6 7 8 9 12 13 14 15 16 17 18],...
    'BDF','off','InitialStep',1e-10,'Stats','off',...
    'AbsTol',1e-12,'RelTol',1e-10,'vectorized','off');
options1 = odeset ('MaxStep',0.00001,'NonNegative',...
    [1 2 3 4 5 6 7 8 9 12 13 14 15 16 17 18],...
    'BDF','off','InitialStep',1e-1,'Stats','off',...
    'AbsTol',1e-14,'RelTol',1e-12,'vectorized','off');
options2 = odeset ('NonNegative',...
    [1 2 3 4 5 6 7 8 9 12 13 14 15 16 17 18],...
    'BDF','on','InitialStep',1e-10,'Stats','off',...
    'AbsTol',1e-12,'RelTol',1e-10);
tic;
% Initial species concentrations
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3 T T_gas 0 0 0 0 0 0 0];
strain = i_strain;
strain_pulse = 0;
tspan = 10000.;
RR_Count = 0;
sol = ode15s(@ammonia4,[0 tspan],s0,options2);
ttt = sol.x';
sss = sol.y';
s0 = sol.y(:,end)';
tspan_cum = tspan;
RTime = toc;
NH3_MF = sss(:,9)./sum(sss(:,7:9),2);
N2_Conv = 1-(1 - NH3_MF)./(1 + NH3_MF);
N2_TOF = c_N2*Q_in*N2_Conv/SDTOT/abyv;
Filename = sprintf('ammonia_strain_Ru_2K_%s_SS.mat',...
                   datestr(now,'mm-dd-yyyy_HH-MM'));
%save(Filename,'-v7.3','ttt','tspan','tspan_cum','s0','period',...
%     'NH3_MF','N2_Conv','N2_TOF','RTime','T_in','SDTOT','c_N2',...
%     'abyv','RATIO_S')
save(Filename,'-v7.3')
T_pulse = T_orig;
%load('ammonia_strain_Ru_SA_Base_200Hz_2steps_02-23-2021_15-56_2.mat','s0','ttt')
%tspan_cum = max(ttt);
%clearvars ttt
strain_pulse = 0;
strain = 0.04;
tspan2 = period/2.;
t_end = tspan2;
sec = 5;
segments = 5;
sol_length = floor(1000*sec/period);
ttt = nan([sol_length 1]);
sss = nan([sol_length 18]);
stst = nan([sol_length 1]);
RR_All = nan(22,3,1500);
tt = nan(1500,1);
rindex = 1;
for count1=1:segments
	RR_Count = 0;
    for count=1:floor(sec/period)
        tspan_cum = tspan_cum + tspan2;
        sol = ode15s(@ammonia4,[0 tspan2],s0,options0);
        t1 = sol.x'+tspan_cum-tspan2;
        s1 = sol.y';
        t_end = max(sol.x');
        rlength  = length(t1);
        ttt(rindex:rindex+rlength-1, :) = t1;
        sss(rindex:rindex+rlength-1, :) = s1;
		st = ones([rlength 1])*strain;
        stst(rindex:rindex+rlength-1, :) = st;
        rindex = rindex + rlength;
        s0 = sol.y(:,end)';
        strain = -strain;
		
        tspan_cum = tspan_cum + tspan2;
        sol = ode15s(@ammonia4,[0 tspan2],s0,options0);
        t2 = sol.x'+tspan_cum-tspan2;
        s2 = sol.y';
        t_end = max(sol.x');
        rlength  = length(t2);
        ttt(rindex:rindex+rlength-1, :) = t2;
        sss(rindex:rindex+rlength-1, :) = s2;
        st = ones([rlength 1])*strain;
        stst(rindex:rindex+rlength-1, :) = st;
        rindex = rindex + rlength;
        s0 = sol.y(:,end)';
        strain = -strain;
		
    end
    RTime = toc;
    end_mat = find(isnan(ttt),1);
    if isempty(end_mat)
        sol_length = floor(length(ttt)*1.2);
    else
        ttt = ttt(1:end_mat-1,:);
        sss = sss(1:end_mat-1,:);
        stst = stst(1:end_mat-1,:);
        sol_length = floor(length(ttt)*1.05);
    end
    NH3_MF = sss(:,9)./sum(sss(:,7:9),2);
    N2_Conv = 1-(1 - NH3_MF)./(1 + NH3_MF);
    N2_TOF = c_N2*Q_in*N2_Conv/SDTOT/abyv;

%	sol_len = length(ttt);
%	RRR = zeros(22,3,sol_len);
%	st = zeros(sol_len, 1);
%	kf = zeros(22, 1, sol_len);
%	kb = zeros(22, 1, sol_len);
%	Ea_hist = zeros(22, 1, sol_len);
%	for ii=1:sol_len
%		strain = nan;
%		[dxdy RRR(:,:,ii) st(ii) kf(:,ii) kb(:,ii) Ea_hist(:,ii)] = ammonia3(ttt(ii),sss(ii,:)');
%	end
    Filename = sprintf('ammonia_strain_Ru_2K_1P_%s_%s.mat',...
        datestr(now,'mm-dd-yyyy_HH-MM'),num2str(count1));
    
%    save(Filename,'-v7.3','ttt','tspan2','tspan_cum','s0','period',...
%        'NH3_MF','N2_Conv','N2_TOF','RTime','T_in','SDTOT','c_N2',...
%	    'abyv','RATIO_S')
    save(Filename,'-v7.3')
    ttt = nan([sol_length 1]);
    sss = nan([sol_length 18]);
    stst = nan([sol_length 1]);
    rindex = 1;
end
tr{1} = sol.x';
sr{1} = sol.y';
NH3_MF = sr{1}(end,9)/sum(sr{1}(end,7:9));
NH3_Conv = (1 - NH3_MF)/(1 + NH3_MF);
Conv = NH3_Conv*100;
NH3_TOF = c_N2*Q_in*(1-NH3_Conv)/SDTOT/abyv;
    fprintf('\n---------------------------------------------\n')
    fprintf('Computation time = %.2f [sec]\n',RTime)
    fprintf('N2 Conversion           = %6.2f [%%]\n',100-NH3_Conv*100)
    fprintf('N2 TOF                  = %6.2e [sec^-1]\n',NH3_TOF)
    fprintf('Temperature (Feed)      = %7.2f [K]\n',T_orig)
    fprintf('Frequency               = %6d [Hz]\n',1/period)
    fprintf('Lattice Strain          = %5.2f %%\n',i_strain*100)
    fprintf('Fraction Step Sites     = %5.2f %%\n',RATIO_S*100)
    fs = sr{1}(end,7:9)./sum(sr{1}(end,7:9),2);
    fprintf('---------------------------------------------\n')
    fprintf('         Gas Species\n')
    fprintf('------------------------------\n')
    fprintf('   N2         H2         NH3\n')
    fprintf('--------   --------   --------\n')
    fprintf('%8.6f   %8.6f   %8.6f\n',fs)
    ss_T = [sr{1}(end,(1:6)) ((SDEN_T*abyv)-sum(sr{1}(end,1:6)))]/...
           (SDEN_T*abyv);
    fprintf('--------------------------------------------------------------------------\n')
    fprintf('                         Terrace Surface Species                          \n')
    fprintf('--------------------------------------------------------------------------\n')
    fprintf(' N2(T)       N(T)       H(T)      NH3(T)     NH2(T)      NH(T)      *(T)\n')
    fprintf('--------   --------   --------   --------   --------   --------   --------\n')
    fprintf('%8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f\n',ss_T)
    fprintf('--------------------------------------------------------------------------\n')
    ss_S = [sr{1}(end,(12:18)) ((SDEN_S*abyv)-sum(sr{1}(end,12:18)))]/(SDEN_S*abyv);
    fprintf('-------------------------------------------------------------------------------------\n')
    fprintf('                                Step Surface Species                                 \n')
    fprintf('-------------------------------------------------------------------------------------\n')
    fprintf(' N2(S1)      N(S)       H(S)      NH3(S)     NH2(S)      NH(S)      N(S3)      *(S)\n')
    fprintf('--------   --------   --------   --------   --------   --------   --------   --------\n')
    fprintf('%8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f\n',ss_S)
    fprintf('-------------------------------------------------------------------------------------\n')
if 0
    figure(1)
    hold on
    plot(tr{1},sr{1}(:,7)./sum(sr{1}(:,7:9),2),'b')
    plot(tr{1},sr{1}(:,8)./sum(sr{1}(:,7:9),2),'r')
    plot(tr{1},sr{1}(:,9)./sum(sr{1}(:,7:9),2),'g')
    hold off
    %xlim([0 tspan+tspan2])
    %xlim([-5 tspan2])
    %ylim([0 1])
    xlabel('Time [sec]')
    ylabel('Mole fraction [gas]')
    legend('N_2','H_2','NH_3')


    figure(2)
    hold on
    plot(tr{1},sr{1}(:,1) ./(SDEN_T*abyv),'-b')
    plot(tr{1},sr{1}(:,2) ./(SDEN_T*abyv),'-r')
    plot(tr{1},sr{1}(:,3) ./(SDEN_T*abyv),'-c')
    plot(tr{1},sr{1}(:,4) ./(SDEN_T*abyv),'-','Color',[0 .45 .74])
    plot(tr{1},sr{1}(:,5) ./(SDEN_T*abyv),'-k')
    plot(tr{1},sr{1}(:,6) ./(SDEN_T*abyv),'-g')
    plot(tr{1},((SDEN_T*abyv)-sum(sr{1}(:,1:6),2))./(SDEN_T*abyv),'-m')
    hold off
    title('Terrace Sites')
    xlim([0 tspan+tspan2])
    % xlim([-5 tspan2])
    ylim([0 1])
    xlabel('Time [sec]')
    ylabel('Surface coverage')
    legend('N_{2*}','N_*','H_*','NH_{3*}','NH_{2*}','NH','\theta_*')

    figure(3)
    hold on
    plot(tr{1},sr{1}(:,12) ./(SDEN_S*abyv),'-b')
    plot(tr{1},sr{1}(:,13) ./(SDEN_S*abyv),'-r')
    plot(tr{1},sr{1}(:,14) ./(SDEN_S*abyv),'-c')
    plot(tr{1},sr{1}(:,15) ./(SDEN_S*abyv),'-','Color',[0 .45 .74])
    plot(tr{1},sr{1}(:,16) ./(SDEN_S*abyv),'-k')
    plot(tr{1},sr{1}(:,17) ./(SDEN_S*abyv),'-g')
    plot(tr{1},sr{1}(:,18) ./(SDEN_S*abyv),'-')
    plot(tr{1},((SDEN_S*abyv)-sum(sr{1}(:,12:18),2))./(SDEN_S*abyv),'-m')
    hold off
    title('Step Sites')
    xlim([0 tspan+tspan2])
    % xlim([-5 tspan2])
    ylim([0 1])
    xlabel('Time [sec]')
    ylabel('Surface coverage')
    legend('N_{2*}','N_*','H_*','NH_{3*}','NH_{2*}','NH','N(S3)','\theta_*')
    hold off

    figure(4)
    TOF = RR ./abyv ./(SDTOT);
    BG = barh(abs(TOF));
    set(BG(3),'DisplayName','Net','FaceColor','r');
    set(BG(2),'DisplayName','Reverse','FaceColor','y');
    set(BG(1),'DisplayName','Forward','FaceColor','b');
    set(gca,'Xscale','log', 'XMinorTick', 'off')
    set(gca,'YTick',1:22)
    set(gca,'yticklabel',{'N_2(T) \leftrightarrow N_2 + *(T)';...
    '2N(T) \leftrightarrow N_2(T) + *(T)';...
    '2H(T) \leftrightarrow H_2 + 2*(T)';...
    'NH(T) + *(T) \leftrightarrow N(T) + H(T)';...
    'NH_2(T) + *(T) \leftrightarrow NH(T) + H(T)';...
    'NH_3(T) + *(T) \leftrightarrow NH_2(T) + H(T)';...
    'NH_3 + *(T) \leftrightarrow NH_3(T)';...
    'N_2(S) \leftrightarrow N_2 + *(S)';...
    '2N(S) \leftrightarrow N_2(S) + *(S)';...
    '2H(S) \leftrightarrow H_2 + 2*(S)';...
    'NH(S) + *(S) \leftrightarrow N(S) + H(S)';...
    'NH_2(S) + *(S) \leftrightarrow NH(S) + H(S)';...
    'NH_3(S) + *(S) \leftrightarrow NH_2(S) + H(S)';...
    'NH_3 + *(S) \leftrightarrow NH_3(S)';...
    'N(T) + *(S) \leftrightarrow N(S) + *(T)';...
    'H(T) + *(S) \leftrightarrow H(S) + *(T)';...
    'NH_3(T) + *(S) \leftrightarrow NH_3(S) + *(T)';...
    'NH_2(T) + *(S) \leftrightarrow NH_2(S) + *(T)';...
    'NH(T) + *(S) \leftrightarrow NH(S) + *(T)';...
    'N_2(S) + *(T) \leftrightarrow N(S) + N(T)';...
    'N_2(S) + *(S) \leftrightarrow N(S) + N(S3)';...
    'N(T) + *(S) \leftrightarrow N(S3) + *(T)';...
    })
    xlabel('Turnover Frequency (TOF) [s^{-1}]')
    ylabel('Reaction Step')
    title({'Ammonia Decomposition', ['Forward and Reverse Reaction Rates at ' ...
            num2str(sr{1}(end,11)) ' [K] on ' Q_name],['V_{Reactor} = ' ...
            num2str(V) ' cm^3     Q_{Feed} = ' num2str(Q_in) ...
            ' cm^3/s     \tau_{Reactor} = ' num2str(V/Q_in) ' seconds']})
    legend('Forward', 'Reverse', 'Net', 'Location', 'best')
    set(gcf, 'position', get(gcf, 'position').*[1 0.64 1.09 1.48])

    figure(5)
    PEI = RR(:,1)./(RR(:,1)+RR(:,2));
    plot(PEI,(1:length(PEI)),'o', 'MarkerFacecolor','b')
    hold on
    h=fill([0.45 0.45 0.55 0.55],[0.5 length(PEI)+0.5 length(PEI)+0.5 0.5],'y');
    set(h,'facealpha',.1,'linestyle','none');
    xlim([0,1])
    set(gca,'YTick',1:22)
    set(gca,'yticklabel',{'N_2(T) \leftrightarrow N_2 + *(T)';...
    '2N(T) \leftrightarrow N_2(T) + *(T)';...
    '2H(T) \leftrightarrow H_2 + 2*(T)';...
    'NH(T) + *(T) \leftrightarrow N(T) + H(T)';...
    'NH_2(T) + *(T) \leftrightarrow NH(T) + H(T)';...
    'NH_3(T) + *(T) \leftrightarrow NH_2(T) + H(T)';...
    'NH_3 + *(T) \leftrightarrow NH_3(T)';...
    'N_2(S) \leftrightarrow N_2 + *(S)';...
    '2N(S) \leftrightarrow N_2(S) + *(S)';...
    '2H(S) \leftrightarrow H_2 + 2*(S)';...
    'NH(S) + *(S) \leftrightarrow N(S) + H(S)';...
    'NH_2(S) + *(S) \leftrightarrow NH(S) + H(S)';...
    'NH_3(S) + *(S) \leftrightarrow NH_2(S) + H(S)';...
    'NH_3 + *(S) \leftrightarrow NH_3(S)';...
    'N(T) + *(S) \leftrightarrow N(S) + *(T)';...
    'H(T) + *(S) \leftrightarrow H(S) + *(T)';...
    'NH_3(T) + *(S) \leftrightarrow NH_3(S) + *(T)';...
    'NH_2(T) + *(S) \leftrightarrow NH_2(S) + *(T)';...
    'NH(T) + *(S) \leftrightarrow NH(S) + *(T)';...
    'N_2(S) + *(T) \leftrightarrow N(S) + N(T)';...
    'N_2(S) + *(S) \leftrightarrow N(S) + N(S3)';...
    'N(T) + *(S) \leftrightarrow N(S3) + *(T)';...
    })
    ylim([0.5, 22.5]);
    plot([0.45 0.45],[0.5,length(PEI)+.5],'--b')
    plot([0.55 0.55],[0.5,length(PEI)+.5],'--b')
    ylabel('Reaction Step')
    xlabel('Partial Equilibrium Index')
    title({'Ammonia Synthesis', ['Partial Equilibrium Index at ' ...
            num2str(sr{1}(end,11)) ' [K] on ' Q_name],['V_{Reactor} = ' ...
            num2str(V) ' cm^3     Q_{Feed} = ' num2str(Q_in) ...
            ' cm^3/s     \tau_{Reactor} = ' num2str(V/Q_in) ' seconds']})
    set(gcf, 'position', get(gcf, 'position').*[1 0.64 1.09 1.48])
    hold off
end
end
