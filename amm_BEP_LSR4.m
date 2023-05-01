function [ Ea,A6_LSR,A6_Strain,A7_Strain,CpOR,HORT,GORT,Q ] = amm_BEP_LSR4( T,Stoic,Ea,strain,A6_Cov)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global R_e Q_name Strain_Coef_H Strain_Coef_S A A_LSR

%% Liner Scaling Relationships

%  Qi = Qi,ref + ALPHAi * (Q_target - Q_ref)

% Use Pt as reference for original thermodynamics data
%Q_ref = 102.35;  % N binding energy on reference metal (Pt) [kcal/mol]
% Use Ru as reference for new thermodynamics data
Q_ref = 134.21;  % N binding energy on reference metal (Ru) [kcal/mol]

%Q_target = 102.35; Q_name = 'Pt'; % N binding energy on Pt
%Q_target = 110.00; Q_name = 'Ni'; % N binding energy on Ni
%Q_target = 112.07; Q_name = 'Rh'; % N binding energy on Rh
%Q_target = 115.30; Q_name = 'Co'; % N binding energy on Co
%Q_target = 132.72; Q_name = 'Os'; % N binding energy on Os
Q_target = 134.21; Q_name = 'Ru'; % N binding energy on Ru
%Q_target = 136.75; Q_name = 'Fe'; % N binding energy on Fe
%Q_target = 138.36; Q_name = 'Re'; % N binding energy on Re
%Q_target = 154.18; Q_name = 'Mo'; % N binding energy on Mo
%Q_target = 102.35; Q_name = 'Unk';

% ALPHAi LSR slopes
alpha(1)  = 0.62036;  %  N2  [Terrace]
alpha(2)  = 1;        %  N   [Terrace]
alpha(3)  = 0.17;     %  H   [Terrace]
alpha(4)  = 0.14;     %  NH3 [Terrace]
alpha(5)  = 0.41;     %  NH2 [Terrace]
alpha(6)  = 0.71;     %  NH  [Terrace]
alpha(7)  = 0.62036;  %  N2  [Step]
alpha(8)  = 1.057;    %  N   [Step]
alpha(9)  = 0.18;     %  H   [Step]
alpha(10) = 0.14;     %  NH3 [Step]
alpha(11) = 0.391;    %  NH2 [Step]
alpha(12) = 0.708;    %  NH  [Step]
alpha(13)  = 1.057;   %  N   [Lower Terrace]

% Qi,ref (zero coverage reference binding energy of the species) [kcal/mol]
Qi_ref(1)  = -2.0779;  %  N2  [Terrace]
Qi_ref(2)  = Q_ref;    %  N   [Terrace]
Qi_ref(3)  = 57.4245;  %  H   [Terrace]
Qi_ref(4)  = 12.2999;  %  NH3 [Terrace]
Qi_ref(5)  = 45.8833;  %  NH2 [Terrace]
Qi_ref(6)  = 82.5372;  %  NH  [Terrace]
Qi_ref(7)  = 9.451;    %  N2  [Step]
Qi_ref(8)  = 106.224;  %  N   [Upper Step]
Qi_ref(9)  = 58.0824;  %  H   [Step]
Qi_ref(10) = 22.6759;  %  NH3 [Step]
Qi_ref(11) = 63.9298;  %  NH2 [Step]
Qi_ref(12) = 91.8554;  %  NH  [Step]
Qi_ref(13)  = 106.224; %  N   [Lower Step]

Q = Qi_ref + alpha * (Q_target - Q_ref);
A6_LSR = (alpha * (Q_target - Q_ref))'/R_e;

%% Catalyst surface strain

A6_Strain = Strain_Coef_H*[strain; 1];
A7_Strain = Strain_Coef_S*[strain^2; strain; 1];
%
%% Bronsted-Evans-Polanyi Relationships for activation barriers from Hrxn

%  (Ea)=m(deltaHrxn)+b

%  m coefficients
m(1) = 0.514;  %N2 dissociation (Terrace)
m(2) = 0.581;  %NH dehydrogenation (Terrace)
m(3) = 0.725;  %NH2 dehydrogenation (Terrace)
m(4) = 0.608;  %NH3 dehydrogenation (Terrace)
m(5) = 0.855;  %N2 dissociation (Step)
m(6) = 0.809;  %NH dehydrogenation (Step)
m(7) = 0.553;  %NH2 dehydrogenation (Step)
m(8) = 0.470;  %NH3 dehydrogenation (Step)
m(9) = 0.183;  %N2 dissociation (Step-Terrace)
m(10)= 0.346;  %N2 dissociation (U Step-L Step)

%  b constant
b(1) = 48.6;   %N2 dissociation (Terrace)
b(2) = 28.0;   %NH dehydrogenation (Terrace)
b(3) = 25.5;   %NH2 dehydrogenation (Terrace)
b(4) = 27.5;   %NH3 dehydrogenation (Terrace)
b(5) = 40.6;   %N2 dissociation (Step)
b(6) = 26.5;   %NH dehydrogenation (Step)
b(7) = 27.7;   %NH2 dehydrogenation (Step)
b(8) = 22.3;   %NH3 dehydrogenation (Step)
b(9) = 18.0;   %N2 dissociation (Step-Terrace)
b(10)= 20.1;   %N2 dissociation (Step-Terrace)

[CpOR,HORT,GORT] = amm_thermo4(T,A6_LSR,A6_Cov,A6_Strain,A7_Strain);
HRXN = HORT * Stoic'*T*R_e;
Ea(2) = m(1) * HRXN(2) + b(1);
Ea(4) = m(4) * HRXN(4) + b(4);
Ea(5) = m(3) * HRXN(5) + b(3);
Ea(6) = m(2) * HRXN(6) + b(2);
Ea(9) = m(5) * HRXN(9) + b(5);
Ea(11) = m(8) * HRXN(11) + b(8);
Ea(12) = m(7) * HRXN(12) + b(7);
Ea(13) = m(6) * HRXN(13) + b(6);
Ea(15) = (strain*100)*(0.46921131) + 20.20056000;  % N* Diffusion
Ea(16) = (strain*100)*(0.10199344) +  9.27356278;  % H* Diffusion
Ea(17) = (strain*100)*(0.06629750) + 13.18263333;  % NH3* Diffusion
Ea(18) = (strain*100)*(0.14124250) +  5.34992000;  % NH2* Diffusion
Ea(19) = (strain*100)*(0.23348250) + 15.03512000;  % NH* Diffusion
Ea(20) = m(9) * HRXN(20) + b(9);
Ea(21) = m(10) * HRXN(21) + b(10);
Ea(22) = (strain*100)*(0.46921131) + 20.20056000;  % N(S3) Diffusion
%% Pre-exponential factor strain scaling

A_New = A_LSR*[strain^2; strain; 1];
A(1:8) = A_New(1:8);
A(14:15) = A_New(9:10);

end