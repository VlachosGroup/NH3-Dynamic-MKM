function [ kf,kb,HORT,CpOR ] = amm_kinetics4( T,T_gas,s )
% Ammonia Synthesis Microkinetic Model
%  Kinetic rate constants
%
% Input:
%   s      Species concentrations
%   T      Temperature, current (K)
%   T_orig Temperature, initial (k)
%
% Output:
%   kf  Forward reaction rate constants
%   kb  Backward reaction rate constants
%
% Species list:
%   s(1)    N2(T)
%   s(2)    N(T)
%   s(3)    H(T)
%   s(4)    NH3(T)
%   s(5)    NH2(T)
%   s(6)    NH(T)
%   s(7)    N2
%   s(8)    H2
%   s(9)    NH3
%   s(10)   Catalyst surface T
%   s(11)   Gas T
%   s(12)   N2(S)
%   s(13)   N(S)
%   s(14)   H(S)
%   s(15)   NH3(S)
%   s(16)   NH2(S)
%   s(17)   NH(S)
%   s(18)   N(S3)
%%
% Kinetic rate constants
global abyv beta Stoic_gas Stoic MWON A Stick MW_N2 MW_H2...
       MW_NH3 R_e R_k R Ea SDTOT strain SDEN_T SDEN_S
T_ref = 1;                      % Reference Temperature [K]
T = T(end);
T_gas = T_gas(end);

%% Lateral Interactions
Effects = zeros(6);
%                      N2*      N*          H*         NH3*        NH2*        NH*         N(S3)
Effects(1:7,1:7) = [-6.00     0.00       0.00  	      0.00        0.00        0.00        0.00;...
                     0.00   -47.0179   -17.7545  	-25.1631    -20.7620    -48.7823    -47.0179;...
                     0.00   -17.7545	-6.7043 	 -9.5019 	 -7.8400    -18.4208    -17.7545;...
                     0.00   -25.1631	-9.5019     -13.4668    -11.1115    -26.1074    -25.1631;...
                     0.00   -20.7620	-7.8400     -11.1115	-19.8316    -21.5412    -20.7620;...
                     0.00   -48.7823   -18.4208     -26.1074   	-21.5412    -50.6129    -48.7823;...
                     0.00   -47.0179   -17.7545     -25.1631    -20.7620    -48.7823    -47.0179];

start_cov = 0;
A6_Cov = zeros(12,1);
A6_Cov(1:6) = (Effects(1:6,1:6)/2*(max(s(1:6)/(SDEN_T*abyv)-start_cov,0)/(1.0-start_cov)))/R_e;
A6_Cov(7:13) = (Effects/2*(max(s(12:18)/(SDEN_S*abyv)-start_cov,0)/(1.0-start_cov)))/R_e;
%%
[Ea,~,~,~,CpOR,HORT,GORT,~] = amm_BEP_LSR4(T,Stoic,Ea,strain,A6_Cov);
kf=zeros(19,1);
%Terrace site reactions
kf(1) = 1*(Stick(1)/(1-MWON*Stick(1)/2))/SDTOT*((T_gas/T_ref)^beta(1))* ...
        sqrt(R_k*T_gas/(2*pi*MW_N2)) * exp(-Ea(1)/(R_e*T_gas));               % N2   +  *(T) <--> N2(T)
kf(2) = 1*A(1)*((T/T_ref)^beta(2))/abyv * exp(-Ea(2)/(R_e*T));                % N2(T)  +  *(T) <--> 2N(T)
kf(3) = 1*(Stick(2)/(1-MWON*Stick(2)/2))/(abyv*SDTOT^2)*...
        ((T_gas/T_ref)^beta(3))*sqrt(R_k*T_gas/(2*pi*MW_H2))*...
        exp(-Ea(3)/(R_e*T_gas));                                              % H2   + 2(T) <--> 2H(T)
kf(4) = 1*A(2)*((T/T_ref)^beta(4))/abyv * exp(-Ea(4)/(R_e*T));                % NH3(T) +  *(T) <--> NH2(T) + H(T)
kf(5) = 1*A(3)*((T/T_ref)^beta(5))/abyv * exp(-Ea(5)/(R_e*T));                % NH2(T) +  *(T) <--> NH(T)  + H(T)
kf(6) = 1*A(4)*((T/T_ref)^beta(6))/abyv * exp(-Ea(6)/(R_e*T));                % NH(T)  +  *(T) <--> N(T)   + H(T)
kf(7) = 1*(Stick(3)/(1-MWON*Stick(3)/2))/SDTOT*((T_gas/T_ref)^beta(7))* ...
        sqrt(R_k*T_gas/(2*pi*MW_NH3)) * exp(-Ea(7)/(R_e*T_gas));              % NH3  +  *(T) <--> NH3(T)

% Step site reactions
kf(8) = 1*(Stick(4)/(1-MWON*Stick(4)/2))/SDTOT*((T_gas/T_ref)^beta(8))* ...
        sqrt(R_k*T_gas/(2*pi*MW_N2)) * exp(-Ea(8)/(R_e*T_gas));               % N2 + *(S) <--> N2(S1)
kf(9) = 0*A(5)*((T/T_ref)^beta(9))/abyv * exp(-Ea(9)/(R_e*T));                % N2(S) + *(S) <--> 2N(S)
kf(10) = 1*(Stick(5)/(1-MWON*Stick(5)/2))/(abyv*SDTOT^2)*...
         ((T_gas/T_ref)^beta(10))*sqrt(R_k*T_gas/(2*pi*MW_H2))*...
         exp(-Ea(10)/(R_e*T_gas));                                            % H2 + 2(S) <--> 2H(S)
kf(11) = 1*A(6)*((T/T_ref)^beta(11))/abyv * exp(-Ea(11)/(R_e*T));             % NH3(S) + *(S) <--> NH2(S) + H(S)
kf(12) = 1*A(7)*((T/T_ref)^beta(12))/abyv * exp(-Ea(12)/(R_e*T));             % NH2(S) + *(S) <--> NH(S) + H(S)
kf(13) = 1*A(8)*((T/T_ref)^beta(13))/abyv * exp(-Ea(13)/(R_e*T));             % NH(S) + *(S) <--> N(S) + H(S)
kf(14) = 1*(Stick(6)/(1-MWON*Stick(6)/2))/SDTOT*((T_gas/T_ref)^beta(14))* ...
         sqrt(R_k*T_gas/(2*pi*MW_NH3)) * exp(-Ea(14)/(R_e*T_gas));            % NH3  +  *(S) <--> NH3(S)
% Diffusion reactions
kf(15) = 1*A(9)*((T/T_ref)^beta(15))/abyv * exp(-Ea(15)/(R_e*T));             % N(T) + *(S) <--> N(S) + *(T)
kf(16) = 0*A(10)*((T/T_ref)^beta(16))/abyv * exp(-Ea(16)/(R_e*T));            % H(T) + *(S) <--> H(S) + *(T)
kf(17) = 1*A(11)*((T/T_ref)^beta(17))/abyv * exp(-Ea(17)/(R_e*T));            % NH3(T) + *(S) <--> NH3(S) + *(T)
kf(18) = 1*A(12)*((T/T_ref)^beta(18))/abyv * exp(-Ea(18)/(R_e*T));            % NH2(T) + *(S) <--> NH2(S) + *(T)
kf(19) = 1*A(13)*((T/T_ref)^beta(19))/abyv * exp(-Ea(19)/(R_e*T));            % NH(T) + *(S) <--> NH(S) + *(T)
kf(20) = 0*A(14)*((T/T_ref)^beta(20))/abyv * exp(-Ea(20)/(R_e*T));            % N2(S) + *(T) <--> N(S) + N(T)
kf(21) = 1*A(15)*((T/T_ref)^beta(21))/abyv * exp(-Ea(21)/(R_e*T));            % N2(S) + *(S) <--> N(S) + N(S3)
kf(22) = 1*A(16)*((T/T_ref)^beta(22))/abyv * exp(-Ea(22)/(R_e*T));            % N(T) + *(S) <--> N(S3) + *(T)

GORT_e = GORT * (Stoic)';
Kp = exp(-GORT_e)';
Kc = Kp .* (1/(R*T)).^(sum(Stoic_gas,2));
kb = kf./Kc;
end

