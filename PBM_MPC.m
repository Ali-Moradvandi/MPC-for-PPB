%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) TU Delft All Rights Reserved
% Written by Abbas Alloul, Giorgio Gardella and Ali Moradvandi
% For any correspondence: a.moradvandi@tudelft.nl

%% Introduction of code (purpose)
% PBM is the main function of the mechanstic model propoesed for purple
% bactria growth in an open raceway reactor. This function  should be
% called in Run_PBM.m file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Structure of the code
% State variables
% Parameters
% Algebraic equations
% Differential equations


function [dydt,Qi,Yxs] = PBM (t,y,Input,Input2,fHS,h)

%% State variables
%  Input composition, flow rates, and paddlewheel
u       = Input(1:15);
SE      = Input2(floor(t)+1,1);
Qin     = Input2(floor(t)+1,2); 
Qout    = Input2(floor(t)+1,3);
sw      = Input2(floor(t)+1,4);
MS      = Input2(floor(t)+1,5);


% Soluble state variables
O2      = y(1);             % Dissolved oxygen                      (mgCOD O2/L)
SS      = y(2);             % Readily biodegradble organic matters  (mgCOD/L)
SVFA    = y(3);             % Volatile fatty acids                  (mgCOD/L)
SIC     = y(4);             % Total inorganic carbon                (mmol HCO3-C/L)
SH2     = y(5);             % Soluble hydrogen                      (mgCOD/L)
SIN     = y(6);             % Total inorganic nitrogen              (mg N/L)
SIP     = y(7);             % Total inorganic phosphorous           (mg P/L)
SI      = y(8);             % Inert soulble organic matters         (mgCOD/L)

% Microbial state variables
XPB_ph  = y(9);             % Photoheterotropic grown PPB           (mgCOD/L)
XPB_ch  = y(10);            % Aerobic chemoheterotropic grown PPB   (mgCOD/L)
XPB_an  = y(11);            % Anaerobic chemoheterotropic grown PPB (mgCOD/L)
XAHB    = y(12);            % Aerobic heterotropic bacteria         (mgCOD/L) 
XAN     = y(13);            % Acidogenic fermenters                 (mgCOD/L)

% Particulate state variables
XS      = y(14);            % Slowly biodegradable organic matter   (mgCOD/L)
XI      = y(15);            % Inert particulate organic matter      (mgCOD/L)

% Volume as a state variable for non-constant volume conditions
Vt      = y(16);            % L

%%  Parameters

% Physical parameters
O2_s        = 7.85;                 % mg O2/L
CO2_s       = 0.0127;               % mmol HCO3-C/L
NH3_s       = 0.533186234;          % mgN/L
H2_s        = 1.6;                  % mg H2/L
Kla_O2      = 1;                    % h-1
Kla_CO2     = 0.9317*sw;            % h-1       % sw is for switching between two modes of paddle wheel.
Kla_NH3     = 0.9165*sw;            % h-1
kla_H2      = 1.5086*sw;            % h-1
Pka_CO2     = 6.37;                 % -
Pka_NH4     = 9.25;                 % -
% Pka_H2PO4   = 7.21;                 % -

% Kinetic parameters - maximal specific growth rates

% MS          = 0.28; 

% if t<=850
%     MS = 0.28;                  % metabolic switch
% % elseif 624<t && t<744
% %     MS = 0.5;
% % elseif 744<t && t<864
% %     MS = 0.05;
% % elseif 864<t && t<984
% %     MS = 0.32;
% % elseif 984<t && t<1104
% %     MS = 0.1;
% else
%     MS = 0.48;
% end

kHYD        = 0.0035;               % h-1
k_dec_PB    = 0.0113;               % h-1
k_dec_AHB   = 0.0156;               % h-1
k_dec_AN    = 0.00083;              % h-1
k_m_VFA_ph  = 0.0783;               % h-1
k_m_SS_ph   = 0.0525;               % h-1
k_m_VFA_ch  = 0.0525;               % h-1
k_m_SS_ch   = 0.05;                 % h-1
k_m_an_PB   = 0.0124;               % h-1
k_m_AHB     = 0.0758;               % h-1
k_m_AN      = 0.0238;               % h-1
   
% Half saturation and inhibition constants
KS_VFA_ph   = 20;                   % mg COD/L
KS_SS_ph    = 5;                    % mg COD/L
KS_VFA_ch   = 0.4;                  % mg COD/L
KS_SS_ch    = 0.4;                  % mg COD/L
KS_SS_an    = 5;                    % mg COD/L
KS_SS_AN    = 5;                    % mg COD/L
KS_VFA_AHB  = 5;                    % mg COD/L
KS_SS_AHB   = 5;                    % mg COD/L
KSIN        = 0.0050;               % mgNH3-N/L
KSIP        = 0.0010;               % mgPO4-P/L
KI_O2_PB    = 0.7;                  % mgCOD O2/L 
KS_O2_PB    = 0.05;                 % mgCOD O2/L 
KI_O2_AHB   = 0.05;                 % mgCOD O2/L 
KI_O2_AN    = 0.05;                 % mgCOD O2/L 

% Biomass Yields - mgCOD/mgCOD
YPBph       = 1.00;                 % PPB photoheterotrophy
YPBch       = 0.52;                 % PPB chemoheterotrophy
YPBan       = 0.197;                % PPB anaerobic chemoheterotrophy
YAHB        = 0.67;                 % AHB chemoheterotrophy
YAN         = 0.197;                % AN anaerobic chemoheterotrophy

% Irradiance parameter
KSE         = 3;                    % W/m2
KIE         = 100;                  % W/m2

% Decay parameters
fICDEC      = -1.98412698412703E-04; % mmol HCO3-C/mgCOD 
fINDEC      = 0.058;                 % mgNH3-N/mgCOD 
fIPDEC      = 0.01;                  % mgPO4-P/mgCOD 

% Hydrolysis parameters
fSSXS       = 1.6382408312061000E-01;   % mgCOD/mgCOD
fSAXS       = 1.166839250294608E-01;    % mgCOD/mgCOD
fICXS       = 1.3039707398869100E-06;   % mmolHCO3-C/mgCOD
fH2XS       = 8.4424680871970500E-02;   % mgCOD/mgCOD
fINXS       = 1.1622464400041300E-02;   % mgNH3-N/mgCOD
fIPXS       = 2.0754400714359400E-03;   % mgPO4-P/mgCOD
fSIXS       = 1.5182085409209700E-01;   % mgCOD/mgCOD
fXIXS       = 4.3309113162071400E-01;   % mgCOD/mgCOD

% Uptake stoichiometry - mmolHCO3-C/mgCOD
fICPHVFA    = 3.897702e-3;          % PPB photoheterotrophy VFA
fICPHSS     = 3.897702e-3;          % PPB photoheterotrophy SS
fICCHVFA    = 0.023728;             % PPB and AHB chemoheterotrophy VFA
fICCHSS     = 0.01400;              % PPB and AHB chemoheterotrophy SS
fICANSS     = -0.02702;             % PPB and AN anaeorbic chemoheterotrophy SS

% AN Chemoheterotrophic VFA formation stochiometry 
fVFAANSS    = 0.7728;               % mgCOD/mgCOD
fH2ANSS     = 0.0304;               % mgCOD/mgCOD

% Nutrient uptake stoichiometry
Y_N_PB     = -0.0860;              % mg N/mg COD
Y_P_PB     = -0.0150;              % mg P/mg COD
Y_N_AHB    = -0.0860;              % mg N/mg COD
Y_P_AHB    = -0.0150;              % mg P/mg COD
Y_N_AN     = -0.0860;              % mg N/mg COD
Y_P_AN     = -0.0150;              % mg P/mg COD

% Oxygen uptake chemoheterotrophy
fO2XPB      = -(1-YPBch)/YPBch;        % mgO2/mgCOD
fO2XAHB     = -(1-YAHB)/YAHB;          % mgO2/mgCOD

%% Algebraic equations
% light, inhibition and rate equations

% light equation: Lambert-Beer s law plus light attenuation
epsilon     = 0.07;         % light extinction coefficien
sigma       = 1.15;         % light factor of absorbance and scattering
SE          = SE*(1 - exp(-epsilon*h*(XPB_ph + XPB_ch + XPB_an+XAHB+XS+XI)/sigma))/...
            (epsilon*h*(XPB_ph + XPB_ch + XPB_an + XAHB + XS + XI)/sigma); % W/m2

PH = 8;

% conversion equations
SCO2        = SIC/(1+(10^(-Pka_CO2))*(10^PH));
SNH3        = SIN/(1+(10^(-PH))/(10^(-Pka_NH4)));


% Inhibition equations
IIN         = SIN/(KSIN + SIN);               % Inorg N limitation
IIP         = SIP/(SIP + KSIP);               % Inorg P limitation
IE_ph       = SE/(SE + KSE);                  % Light inhibition
IE_ch       = 1-(SE/(SE + KIE));              % Light inhibition chemo
IO2_XPB     = 1-(O2/(KI_O2_PB + O2));         % Oxygen inhibithion PPB photoheterotrophy
IO2_XPB2    = O2/(KS_O2_PB + O2);             % Oxygen inhibithion PPB chemoheterotrophy 
IO2_AHB     = O2/(KI_O2_AHB + O2);            % Oxygen affinity AHB
IO2_AN      = 1-(O2/(KI_O2_AN + O2));         % Oxygen inhibithion AN and PPB fermentation
ICs         = SVFA/(SS + SVFA);               % Competitive inhibition PPB on VFA due to SS
ICVFA       = SS/(SS + SVFA);                 % Competitive inhibition PPB on VFA due to SS

% Rate equations

XPB         = XPB_ph + XPB_ch + XPB_an;
rHYD        = kHYD*XS;                                                      % Hydrolysis
r_dec_PB    = k_dec_PB*XPB;                                                 % Decay PPB
r_dec_AHB   = k_dec_AHB*XAHB;                                               % Decay AHB
r_dec_AN    = k_dec_AN*XAN;                                                 % Decay AN
r_VFA_ph    = k_m_VFA_ph*(XPB_ph+XPB_ch*MS + XPB_an*MS)*IE_ph*IO2_XPB...
             *IIN*IIP*ICs*(SVFA/(KS_VFA_ph + SVFA));                        % PPB photoheterotrophy VFA
r_ss_ph     = k_m_SS_ph *(XPB_ph+XPB_ch*MS + XPB_an*MS)*IE_ph*IO2_XPB...
             *IIN*IIP*ICVFA*(SS/(KS_SS_ph + SS));                           % PPB photoheterotrophy SS
r_VFA_ch    = k_m_VFA_ch* (XPB_ph*MS+XPB_ch+ XPB_an*MS)*(IE_ch)*IO2_XPB2...
             *IIN*IIP *ICs*(SVFA/(KS_VFA_ch + SVFA));                       % PPB chemoheterotrophy VFA
r_ss_ch     = k_m_SS_ch * (XPB_ph*MS+XPB_ch+ XPB_an*MS)*(IE_ch)*IO2_XPB2...
             *IIN*IIP*ICVFA*(SS/(KS_SS_ch + SS));                           % PPB chemoheterotrophy SS
r_ss_an_pb  = k_m_an_PB * (XPB_ph*MS+XPB_ch*MS+ XPB_an) * (IE_ch) *IO2_AN...
             *IIN*IIP*(SS/(KS_SS_an + SS));                                 % PPB anaerobic chemoheterotrophy SS
r_VFA_AHB   = k_m_AHB*XAHB*IO2_AHB*IIN*IIP*ICs*(SVFA/(KS_VFA_AHB + SVFA));  % AHB chemoheterotrophy VFA
r_ss_AHB    = k_m_AHB*XAHB*IO2_AHB*IIN*IIP*ICVFA*(SS/(KS_SS_AHB + SS));     % AHB chemoheterotrophy SS
r_ss_an     = k_m_AN*XAN*IO2_AN*IIN*IIP*(SS/(KS_SS_AN + SS));               % AN anaerobic chemoheterotrophy SS


r_tot       = r_VFA_ph + r_ss_ph + r_VFA_ch + r_ss_ch + r_ss_an_pb; 
r_tot_AHB   = r_VFA_AHB + r_ss_AHB;
r_dec_tot   = r_dec_PB + r_dec_AHB + r_dec_AN;

%% Differential equations

dydt     = zeros(16,1);

% Volume
dydt(16) = Qin - Qout;

% O2
dydt(1) = (u(1)*Qin - O2*Qout - dydt(16)*O2)/Vt +(Kla_O2*(O2_s - O2)...
        +fO2XPB*r_ss_ch + fO2XPB*r_VFA_ch+ fO2XAHB*r_ss_AHB...
        + fO2XAHB*r_VFA_AHB)*sw - Kla_O2*O2*(1-sw);
% SS    
dydt(2) = (u(2)*Qin - SS*Qout - dydt(16)*SS)/Vt...
        -r_ss_ph/YPBph -r_ss_ch/YPBch -r_ss_an_pb/YPBan...
        -r_ss_AHB/YAHB - r_ss_an/YPBan + fSSXS * rHYD; 
% SVFA
dydt(3) = (u(3)*Qin - SVFA*Qout - dydt(16)*SVFA)/Vt - r_VFA_ph/YPBph...
        -r_VFA_ch/YPBch - r_VFA_AHB/YAHB + r_ss_an/YAN*fVFAANSS ...
        + r_ss_an_pb/YAN*fVFAANSS + fSAXS*rHYD;
% SIC
dydt(4) = (u(4)*Qin - (SIC)*Qout - dydt(16)*(SIC))/Vt + fICPHVFA * r_VFA_ph...
        + fICPHSS * r_ss_ph + fICCHVFA*r_VFA_ch +fICCHSS*r_ss_ch...
        + fICANSS * r_ss_an_pb+ fICCHVFA*r_VFA_AHB + fICCHSS*r_ss_AHB...
        + fICANSS*r_ss_an + fICDEC*r_dec_tot ...
        + Kla_CO2*(CO2_s-SCO2)+ fICXS*rHYD ;                       
% H2
dydt(5) = (u(5)*Qin - SH2*Qout- dydt(16)*SH2)/Vt...
        + fH2XS * rHYD + kla_H2*(H2_s-SH2) + r_ss_an/YAN*fH2ANSS;
% SIN
dydt(6) =(u(6)*Qin - SIN*Qout - dydt(16)*SIN)/Vt + Y_N_PB*(r_tot)...
        + Y_N_AHB*(r_tot_AHB) + Y_N_AN*r_ss_an + fINXS*rHYD...
        + fINDEC*r_dec_tot +Kla_NH3*(NH3_s-SNH3);             
% SIP
dydt(7) = (u(7)*Qin - SIP*Qout - dydt(16)*SIP)/Vt+Y_P_PB*(r_tot)...
        + Y_P_AHB*(r_tot_AHB) + Y_P_AN * r_ss_an + fIPXS*rHYD...
        + fIPDEC*r_dec_tot;  
% SI
dydt(8) = (u(8)*Qin - SI*Qout-dydt(16)*SI)/Vt + fSIXS*rHYD;
% XPB_ph
dydt(9) = (u(9)*Qin - XPB_ph*fHS*Qout - dydt(16)*XPB_ph)/Vt...
        + r_VFA_ph + r_ss_ph - r_dec_PB/XPB*XPB_ph;
% XPB_ch
dydt(10) = (u(10)*Qin - XPB_ch*fHS*Qout - dydt(16)*XPB_ch)/Vt...
        + r_VFA_ch + r_ss_ch - r_dec_PB/XPB*XPB_ch;
% XPB_an
dydt(11) = (u(11)*Qin - XPB_an*fHS*Qout - dydt(16)*XPB_an)/Vt...
        + r_ss_an_pb - r_dec_PB/XPB*XPB_an;
% XAHB
dydt(12) = (u(12)*Qin - XAHB*fHS*Qout - dydt(16)*XAHB)/Vt ...
        + r_tot_AHB - r_dec_AHB; 
% XAN
dydt(13) = (u(13)*Qin - XAN*fHS*Qout - dydt(16)*XAN)/Vt ...
        + r_ss_an - r_dec_AN;
% XS
dydt(14) = (u(14)*Qin - XS*fHS*Qout - dydt(16)*XS)/Vt - rHYD...
        + r_dec_tot;               
% XI
dydt(15) = (u(15)*Qin - XI*fHS*Qout - dydt(16)*XI)/Vt + fXIXS*rHYD;
 

%% Outputs
Qi      = [Qin, Qout, dydt(16)];
Y_PB    = (r_ss_ph +r_ss_ch+r_VFA_ph +r_VFA_ch)/...
         (r_ss_ph/YPBph + r_VFA_ph/YPBph + r_ss_ch/YPBch + r_VFA_ch/YPBch);
Y_AHB   = (r_ss_AHB + r_VFA_AHB)/(r_ss_AHB/YAHB + r_VFA_AHB/YAHB);
Y_TOT   = (XPB+XAHB+XAN)/(XPB/Y_PB + XAHB/Y_AHB + XAN/YAN);
Yxs     = [Y_PB,Y_AHB,Y_TOT];

end