function [ CpOR,HORT,GORT ] = amm_thermo4( T,A6_LSR,A6_Cov,A6_Strain,A7_Strain)
%%          -------------------------------------------------
%                        NH3  Micro-kinetic model
%                         Vlachos Research Group
%                 Chemical and Biomolecular Egineering
%                         University of Delaware
%
%             Gerhard R Wittreich, P.E.  (February 10, 2017)
%           --------------------------------------------------
%
%      amm_thermo.m   : NASA polynomials provides enthalpy and entropy
%
%      requires: 
%                amm_kinetics.m : Provides forward and reverse rate constants
%
% Input:
%   T     Temperature
%   A6_LSR  LSR correction to NASA polynomial A6 coefficient
%   A6_COV  Coverage correction to NASA polynomial A6 coefficient
%
% Output:
%   HORT   Dimmensionless enthalpy
%   SOR    Dimmensionless entropy
%   GORT   Dimmensionless Gibb's free energy
%
T_c = 500;
T = T(end);
%% pMuTT Recalculated Terrace site thermodynamic data w/ Shizhong DFT calculations Tref = 298.15 K
A_N2s1_h =  [ 3.88555607E+00  2.25579107E-03 -6.74243181E-07 -1.56805094E-10  8.86546530E-14 -7.99940790e+03 -1.58249660e+01];
A_N2s1_l =  [-4.15269022E-01  3.85284631E-02 -1.17383449E-04  1.68785681E-07 -9.23186795E-11 -7.58237629e+03  1.75988553e+00];
A_Ns1_h =   [ 8.81592750E-01  5.93998550E-03 -6.89025148E-06  3.70304607E-09 -7.59376769E-13 -1.32276911e+04 -6.40687724e+00];
A_Ns1_l =   [-6.26549307E-01  3.10362616E-03  5.83736827E-05 -1.81844809E-07  1.56189066E-10 -1.29201483e+04  1.50466373e+00];
A_Hs1_h =   [-1.59446486E+00  1.16296966E-02 -1.25644676E-05  6.41488372E-09 -1.26617068E-12 -6.34177280e+03  6.18405729e+00];
A_Hs1_l =   [ 6.88496466E-01 -1.38512868E-02  8.47127091E-05 -1.51068877E-07  9.17686332E-11 -6.47213000e+03 -2.31465351e+00];
A_NH3s1_h = [ 2.59608053e+00  1.05184631e-02 -7.12704042e-06  3.18299981e-09 -6.48423332e-13 -1.60457863e+04 -1.09600590e+01];
A_NH3s1_l = [ 9.36635813e-01  2.23568480e-02 -4.08447870e-05  4.76221396e-08 -2.29237157e-11 -1.58460969e+04 -3.85531391e+00];
A_NH2s1_h = [ 1.72904798e+00  1.12305142e-02 -1.05353687e-05  5.52410610e-09 -1.16990017e-12 -1.40390955e+04 -1.05499510e+01];
A_NH2s1_l = [-5.52107493e-01  1.28669584e-02  4.03721594e-05 -1.51178033e-07  1.31909234e-10 -1.36074939e+04  8.94767559e-01];
A_NHs1_h =  [ 1.12030379e+00  9.97260213e-03 -1.08784221e-05  5.96143775e-09 -1.26605718e-12 -1.68745942e+04 -8.56185727e+00];
A_NHs1_l =  [-6.18882504e-01  1.00054202e-03  9.68705805e-05 -2.70604612e-07  2.20447297e-10 -1.64373993e+04  1.32322291e+00];

%% pMuTT Recalculated Step site thermodynamic data w/ Shizhong DFT calculations Tref = 298.15 K
A_N2s2_h =  [ 1.97820282e+00  1.07131158e-02 -1.20166339e-05  6.31028436e-09 -1.27266149e-12 -1.03805148e+04 -1.16443181e+01];
A_N2s2_l =  [-2.81331763e+00  4.24076141e-02 -8.66566041e-05  7.36959255e-08 -1.44374379e-11 -9.80718890e+03  9.01381543e+00];
A_Ns2_h =   [ 8.41886313e-01  6.04605571e-03 -7.00922981e-06  3.76545786e-09 -7.71945940e-13 -1.38038997e+04 -6.24521075e+00];
A_Ns3_l =   [-9.90256114E-01  9.89041505E-03  2.66141929E-05 -1.21466365E-07  1.14826904E-10 -7.99837762E+03  2.79729570E+00];
A_Ns3_h =   [ 1.10948357E+00  5.32401108E-03 -6.19354861E-06  3.33534788E-09 -6.84984710E-13 -8.33853519E+03 -7.26279199E+00];
A_Ns2_l =   [-5.67223056e-01  2.06423227e-03  6.29862277e-05 -1.90213926e-07  1.61692659e-10 -1.35025704e+04  1.29729614e+00];
A_Hs2_h =   [-9.77206800e-01  8.99327950e-03 -8.89490036e-06  4.23588302e-09 -7.90961070e-13 -5.95165410e+03  4.39039306e+00];
A_Hs2_l =   [-6.19376989e-01  1.22983220e-02 -4.41535494e-05  9.62968227e-08 -7.61423868e-11 -6.04209472e+03  2.26292442e+00];
A_NH3s2_h = [ 2.62924993e+00  1.03616497e-02 -6.89457817e-06  3.04147991e-09 -6.17352171e-13 -1.89502630e+04 -1.11844959e+01];
A_NH3s2_l = [ 8.56500778e-01  2.34498342e-02 -4.53414530e-05  5.50424742e-08 -2.72785395e-11 -1.87438414e+04 -3.65591455e+00];
A_NH2s2_h = [ 1.82458452e+00  1.10967319e-02 -1.04980610e-05  5.54164275e-09 -1.17777445e-12 -2.05902631e+04 -1.06632044e+01];
A_NH2s2_l = [-2.03257483e-01  1.05402360e-02  4.82281444e-05 -1.64439376e-07  1.40688407e-10 -2.01844156e+04 -2.57605468e-01];
A_NHs2_h =  [ 1.31532728e+00  9.44907643e-03 -1.02762378e-05  5.64002517e-09 -1.20073535e-12 -1.57710978e+04 -9.38947309e+00];
A_NHs2_l =  [-8.23906053e-01  4.78422346e-03  8.05290632e-05 -2.41782079e-07  2.01944194e-10 -1.53056142e+04  2.02190235e+00];

%% NIST gas phase thermodynamic data Tref = 298.15 K
A_H2_h =  [ 3.87477067E+00 -1.37112242E-03  1.67715317E-06 -6.29819288E-10  8.25486826E-14 -1.12008569E+03 -6.05293502E+00];
A_H2_l =  [ 2.52847231E+00  7.09313422E-03 -1.91040770E-05  2.27694673E-08 -1.00153522E-11 -9.40599211E+02 -1.36077347E-01];
A_N2_h =  [ 2.95777821E+00  1.13094739E-03  3.65148606E-08 -2.56570523E-10  6.33874587E-14 -9.02652143E+02  5.93131992E+00];
A_N2_l =  [ 3.47021080E+00  5.16439891E-04 -2.88494917E-06  6.06421706E-09 -3.32545322E-12 -1.04251991E+03  3.20064728E+00];
A_NH3_h = [ 2.78395662E+00  4.95809145E-03 -6.99897857E-07 -3.64861853E-10  1.18192265E-13 -6.58708192E+03  5.86245616E+00];
A_NH3_l = [ 4.38591560E+00 -6.45610533E-03  3.06335822E-05 -3.93962796E-08  1.86235890E-11 -6.77288321E+03 -9.65971735E-01];
A_v_h = [0 0 0 0 0 0 0];
A_v_l = [0 0 0 0 0 0 0];

%%
if T>791.5
    A_H2 = A_H2_h;
    A_N2 = A_N2_h;
else
    A_H2 = A_H2_l;
    A_N2 = A_N2_l;
end
if T>592.4
    A_NH3 = A_NH3_h;
else
    A_NH3 = A_NH3_l;
end
A_h = [A_N2s1_h;A_Ns1_h;A_Hs1_h;A_NH3s1_h;A_NH2s1_h;A_NHs1_h;A_N2;A_H2;A_NH3;A_v_h;...
       A_N2s2_h;A_Ns2_h;A_Hs2_h;A_NH3s2_h;A_NH2s2_h;A_NHs2_h;A_v_h;A_Ns3_h];
A_l = [A_N2s1_l;A_Ns1_l;A_Hs1_l;A_NH3s1_l;A_NH2s1_l;A_NHs1_l;A_N2;A_H2;A_NH3;A_v_l;...
       A_N2s2_l;A_Ns2_l;A_Hs2_l;A_NH3s2_l;A_NH2s2_l;A_NHs2_l;A_v_l;A_Ns3_l];

if T>T_c
    A = A_h;
else
    A = A_l;
end
A6_Correction = A6_LSR + A6_Cov - A6_Strain;
A(1:6,6) =  A(1:6,6) - A6_Correction(1:6);
A(11:16,6) =  A(11:16,6) - A6_Correction(7:12);
A(18,6) = A(18,6) - A6_Correction(13);
A(1:6,7) =  A(1:6,7) + A7_Strain(1:6);
A(11:16,7) =  A(11:16,7) + A7_Strain(7:12);
A(18,7) = A(18,7) + A7_Strain(13);
T_Cp = [1 T T^2 T^3 T^4];
T_H  = [1 T/2 T^2/3 T^3/4 T^4/5 1/T];
T_S  = [log(T) T T^2/2 T^3/3 T^4/4 1];
CpOR = T_Cp*A(:,1:5)';
HORT = T_H*A(:,1:6)';
SOR  = T_S*A(:,[1:5 7])';
GORT = HORT - SOR;
end

