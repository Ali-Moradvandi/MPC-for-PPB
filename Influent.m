%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) TU Delft All Rights Reserved
% Written by Abbas Alloul, Giorgio Gardella and Ali Moradvandi
% For any correspondence: a.moradvandi@tudelft.nl

%% Introduction of code (purpose)
% it is for construction of the influent.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Influent] = Influent 

%% State variable   = value         % Unit
O2                  = 0;            % mg O2/L
SS                  = 0;            % mgCODSS/L
SVFA                = 3000;         % mgCODHAC/L
SIC                 = 3.57;         % mmol HCO3-C/L
SH2                 = 0;            % mgCOD H2 /L
SIN                 = 314;          % mgNH4-N/L
SIP                 = 182;          % mgPO4-P/L
SI                  = 0;            % mgCOD/L
XPB_ph              = 0;            % mg CODx/L 1.75gCOD/gTSS = 20mgTSS/L
XPB_ch              = 0;            % mg CODx/L
XPB_an              = 0;            % mg CODx/L
XAHB                = 0;            % mgCOD/L
XAN                 = 0;            % mgCOD/L
XS                  = 0;            % mgCOD/L
XI                  = 0;            % mgCOD/L


Influent = [O2,SS,SVFA,SIC,SH2,SIN,SIP,SI,XPB_ph,XPB_ch,...
    XPB_an,XAHB,XAN,XS,XI];

end 

