function [ ds RR strain kf kb Ea] = ammonia4( t,s )
%%          -------------------------------------------------
%                        NH3  Micro-kinetic model
%                         Vlachos Research Group
%                 Chemical and Biomolecular Egineering
%                         University of Delaware
%
%             Gerhard R Wittreich, P.E.  (February 10, 2017)
%           --------------------------------------------------
%
%
%  ammonia.m      : Reaction ODE's and dynamical equations
%
%      requires: 
%                amm_kinetics.m : Provides forward and reverse rate constants
%
% Input:
%   s   Species concentrations
%   T   Temperature (K)
%
% Output:
%   ds  Reaction rate for each species
%
% Species list:
%   [s(1),ds(1)]    N2(T)
%   [s(2),ds(2)]    N(T)
%   [s(3),ds(3)]    H(T)
%   [s(4),ds(4)]    NH3(T)
%   [s(5),ds(5)]    NH2(T)
%   [s(6),ds(6)]    NH(T)
%   [s(7),ds(7)]    N2
%   [s(8),ds(8)]    H2
%   [s(9),ds(9)]    NH3
%   [vac_T,dvac_T]   * (Vacant site)
%   [s(10),ds(10)]  Catalyst surface temperature
%   [s(11),ds(11)]  Gas temperature
%   [s(12),ds(12)]  N2(S)
%   [s(13),ds(13)]  N(S)
%   [s(14),ds(14)]  H(S)
%   [s(15),ds(15)]  NH3(S)
%   [s(16),ds(16)]  NH2(S)
%   [s(17),ds(17)]  NH(S)
%   [s(18),ds(18)]  N(S3)
%   [vac_S,dvac_s]  *(S) (Vacant site)
%
global T_gas V Q_in c_N2 c_H2 c_NH3 Isobaric t_check...
       abyv RR strain_pulse strain tspan  SDEN_T SDEN_S period tspan_cum...
	   RR_Count RR_All tt Ea
RR_Count = RR_Count + 1;
T_gas = s(11);
T = s(10);
vac_T = (SDEN_T*abyv) - sum(s(1:6));
vac_S = (SDEN_S*abyv) - sum(s(12:18));
if isnan(strain)
	magnitude = 0.04; % Pulse amplitude
	strain = square(2*pi*t/period)*magnitude;
end
if strain_pulse
    magnitude = 0.04; % Pulse amplitude
    strain = sin(2*pi*t/period)*magnitude;
%    if t <= tspan+period/4.
%        strain = sin(2*pi*(t-tspan)/period)^3*magnitude;
%    else
%        strain = sin(2*pi*(t-tspan)/period)*magnitude;
%    end
end

[kf,kb,~,~]=amm_kinetics4(T,T_gas,s);  % Obtain kinetics rate constants

% Adjust reactor outflow to maintain an isobaric reactor if Isobaric = 1
Q_r = (kb(1)*s(1)     - kf(1)*s(7)*vac_T +...
       kb(3)*s(3)^2   - kf(3)*s(8)*vac_T^2 +...
       kb(7)*s(4)     - kf(7)*s(9)*vac_T +...
       kb(8)*s(12)    - kf(8)*s(7)*vac_S +...
       kb(10)*s(14)^2 - kf(10)*s(8)*vac_S^2 +...
       kb(14)*s(15)   - kf(14)*s(9)*vac_S)/sum(s(7:9))*V;
Q_out = Q_in + Q_r*Isobaric;
%
% Reaction network
%
RR = zeros(22,3);
RR(1,1)  = kf(1)*s(7)*vac_T;
RR(1,2)  = kb(1)*s(1);
RR(1,3)  = RR(1,1) - RR(1,2);
RR(2,1)  = kf(2)*s(1)*vac_T;
RR(2,2)  = kb(2)*s(2)^2;
RR(2,3)  = RR(2,1) - RR(2,2);
RR(3,1)  = kf(3)*s(8)*vac_T^2;
RR(3,2)  = kb(3)*s(3)^2;
RR(3,3)  = RR(3,1) - RR(3,2);
RR(4,1)  = kf(6)*s(6)*vac_T;
RR(4,2)  = kb(6)*s(2)*s(3);
RR(4,3)  = RR(4,1) - RR(4,2);
RR(5,1)  = kf(5)*s(5)*vac_T;
RR(5,2)  = kb(5)*s(6)*s(3);
RR(5,3)  = RR(5,1) - RR(5,2);
RR(6,1)  = kf(4)*s(4)*vac_T;
RR(6,2)  = kb(4)*s(5)*s(3);
RR(6,3)  = RR(6,1) - RR(6,2);
RR(7,1)  = kf(7)*s(9)*vac_T;
RR(7,2)  = kb(7)*s(4);
RR(7,3)  = RR(7,1) - RR(7,2);
RR(8,1)  = kf(8)*s(7)*vac_S;
RR(8,2)  = kb(8)*s(12);
RR(8,3)  = RR(8,1) - RR(8,2);
RR(9,1)  = kf(9)*s(12)*vac_S;
RR(9,2)  = kb(9)*s(13)^2;
RR(9,3)  = RR(9,1) - RR(9,2);
RR(10,1) = kf(10)*s(8)*vac_S^2;
RR(10,2) = kb(10)*s(14)^2;
RR(10,3) = RR(10,1) - RR(10,2);
RR(11,1) = kf(11)*s(15)*vac_S;
RR(11,2) = kb(11)*s(16)*s(14);
RR(11,3) = RR(11,1) - RR(11,2);
RR(12,1) = kf(12)*s(16)*vac_S;
RR(12,2) = kb(12)*s(17)*s(14);
RR(12,3) = RR(12,1) - RR(12,2);
RR(13,1) = kf(13)*s(17)*vac_S;
RR(13,2) = kb(13)*s(13)*s(14);
RR(13,3) = RR(13,1) - RR(13,2);
RR(14,1) = kf(14)*s(9)*vac_S;
RR(14,2) = kb(14)*s(15);
RR(14,3) = RR(14,1) - RR(14,2);
RR(15,1) = kf(15)*s(2)*vac_S;
RR(15,2) = kb(15)*s(13)*vac_T;
RR(15,3) = RR(15,1) - RR(15,2);
% RR(16,1) = kf(16)*s(3)*vac_S;
% RR(16,2) = kb(16)*s(14)*vac_T;
% RR(16,3) = RR(16,1) - RR(16,2);
RR(17,1) = kf(17)*s(4)*vac_S;
RR(17,2) = kb(17)*s(15)*vac_T;
RR(17,3) = RR(17,1) - RR(17,2);
RR(18,1) = kf(18)*s(5)*vac_S;
RR(18,2) = kb(18)*s(16)*vac_T;
RR(18,3) = RR(18,1) - RR(18,2);
RR(19,1) = kf(19)*s(6)*vac_S;
RR(19,2) = kb(19)*s(17)*vac_T;
RR(19,3) = RR(19,1) - RR(19,2);
RR(20,1) = kf(20)*s(12)*vac_T;
RR(20,2) = kb(20)*s(2)*s(13);
RR(20,3) = RR(20,1) - RR(20,2);
RR(21,1) = kf(21)*s(12)*vac_S;
RR(21,2) = kb(21)*s(13)*s(18);
RR(21,3) = RR(21,1) - RR(21,2);
RR(22,1) = kf(22)*s(2)*vac_S;
RR(22,2) = kb(22)*s(18)*vac_T;
RR(22,3) = RR(22,1) - RR(22,2);

ds    = zeros(size(s),'like',s);
ds(1) = kf(1)*s(7)*vac_T         - kb(1)*s(1) + ...
        kb(2)*s(2)^2             - kf(2)*s(1)*vac_T;             % dN2(T)/dt
ds(2) = 2*kf(2)*s(1)*vac_T       - 2*kb(2)*s(2)^2 + ...
        kf(6)*s(6)*vac_T         - kb(6)*s(2)*s(3) +...
        kb(15)*s(13)*vac_T       - kf(15)*s(2)*vac_S +...
        kf(20)*s(12)*vac_T       - kb(20)*s(2)*s(13) +...
        kb(22)*s(18)*vac_T       - kf(22)*s(2)*vac_S;            % dN(T)/dt
ds(3) = 2*kf(3)*s(8)*vac_T^2     - 2*kb(3)*s(3)^2 +...
        kf(4)*s(4)*vac_T         - kb(4)*s(5)*s(3) +...
        kf(5)*s(5)*vac_T         - kb(5)*s(6)*s(3) +...
        kf(6)*s(6)*vac_T         - kb(6)*s(2)*s(3) + ...
        kb(16)*s(14)*vac_T       - kf(16)*s(3)*vac_S;            % dH(T)/dt
ds(4) = kb(4)*s(5)*s(3)          - kf(4)*s(4)*vac_T +...
        kf(7)*s(9)*vac_T         - kb(7)*s(4) + ...
        kb(17)*s(15)*vac_T       - kf(17)*s(4)*vac_S;            % dNH3(T)/dt
ds(5) = kf(4)*s(4)*vac_T         - kb(4)*s(5)*s(3) +...
        kb(5)*s(6)*s(3)          - kf(5)*s(5)*vac_T +...
        kb(18)*s(16)*vac_T       - kf(18)*s(5)*vac_S;            % dNH2(T)/dt
ds(6) = kf(5)*s(5)*vac_T         - kb(5)*s(6)*s(3) +...
        kb(6)*s(2)*s(3)          - kf(6)*s(6)*vac_T +...
        kb(19)*s(17)*vac_T       - kf(19)*s(6)*vac_S;            % dNH(T)/dt
ds(7) = kb(1)*s(1)               - kf(1)*s(7)*vac_T +...
        kb(8)*s(12)              - kf(8)*s(7)*vac_S +...
        Q_in/V*c_N2 - Q_out/V*s(7);                              % dN2/dt
ds(8) = kb(3)*s(3)^2             - kf(3)*s(8)*vac_T^2 +...
        kb(10)*s(14)^2           - kf(10)*s(8)*vac_S^2 +...
        Q_in/V*c_H2 - Q_out/V*s(8);                              % dH2/dt
ds(9) = kb(7)*s(4)               - kf(7)*s(9)*vac_T +...
        kb(14)*s(15)             - kf(14)*s(9)*vac_S +...
        Q_in/V*c_NH3 - Q_out/V*s(9);                             % dNH3/dt
ds(12) = kf(8)*s(7)*vac_S        - kb(8)*s(12) + ...
         kb(9)*s(13)^2           - kf(9)*s(12)*vac_S +...
         kb(20)*s(2)*s(13)       - kf(20)*s(12)*vac_T +...
         kb(21)*s(13)*s(18)      - kf(21)*s(12)*vac_S;           % dN2(S)/dt
ds(13) = 2*kf(9)*s(12)*vac_S     - 2*kb(9)*s(13)^2 + ...
         kf(13)*s(17)*vac_S      - kb(13)*s(13)*s(14) +...
         kf(15)*s(2)*vac_S       - kb(15)*s(13)*vac_T +...
         kf(20)*s(12)*vac_T      - kb(20)*s(2)*s(13) +...
         kf(21)*s(12)*vac_S      - kb(21)*s(13)*s(18);           % dN(S)/dt
ds(14) = 2*kf(10)*s(8)*vac_S^2   - 2*kb(10)*s(14)^2 +...
         kf(11)*s(15)*vac_S      - kb(11)*s(16)*s(14) +...
         kf(12)*s(16)*vac_S      - kb(12)*s(17)*s(14) +...
         kf(13)*s(17)*vac_S      - kb(13)*s(13)*s(14) +...
         kf(16)*s(3)*vac_S       - kb(16)*s(14)*vac_T;           % dH(S)/dt
ds(15) = kb(11)*s(16)*s(14)      - kf(11)*s(15)*vac_S +...
         kf(14)*s(9)*vac_S       - kb(14)*s(15) + ...
         kf(17)*s(4)*vac_S       - kb(17)*s(15)*vac_T;           % dNH3(S)/dt
ds(16) = kf(11)*s(15)*vac_S      - kb(11)*s(16)*s(14) +...
         kb(12)*s(17)*s(14)      - kf(12)*s(16)*vac_S +...
         kf(18)*s(5)*vac_S       - kb(18)*s(16)*vac_T;           % dNH2(S)/dt
ds(17) = kf(12)*s(16)*vac_S      - kb(12)*s(17)*s(14) +...
         kb(13)*s(13)*s(14)      - kf(13)*s(17)*vac_S +...
         kf(19)*s(6)*vac_S       - kb(19)*s(17)*vac_T;           % dNH(S)/dt
ds(18) = kf(21)*s(12)*vac_S      - kb(21)*s(13)*s(18) +...
         kf(22)*s(2)*vac_S       - kb(22)*s(18)*vac_T;           % dN(S3)/dt
%if floor(t*10)>t_check*10 && t>tspan
%    t_check = floor(t*10)/10;
%	fprintf('%4.1f\n',t_check)
%end
    
end
