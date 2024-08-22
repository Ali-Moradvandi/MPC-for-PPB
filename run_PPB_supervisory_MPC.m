%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) TU Delft All Rights Reserved
% Written by Ali Moradvandi
% For any correspondence: moradvandi@gmail.com

%% Introduction of code (purpose)
% This code is the main code for the control system of purple
% phototrophic bacteria (PPB) growth in raceway reactors based on Model
% Predictive Control (MPC).
% The details of the control system configuration can be found
% in the paper archived in https://ssrn.com/abstract=4921489

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model predictive control for a raceway reactor to cultivate PPB

clc 
clear all

% load initial condition from the PBM model for PPB-dominant condition
% PBM should run with paddlewheel off for 50 days with HRT=SRT=4days
load init_pdoff_50d_44.mat  
initial     = y(960,:); 

% for manual adjusting the initial condition you may use:
% Initial conditions % [O2, SS, SVFA, SIC, SH2, SIN, SIP, SI, XPB_ph, XPB_ch, XPB_an, XAHB, XAN, XS, XI, V] 
% initial     = [0, 0, 3000, 3.57, 0, 314,  182 ,0, 80, 10, 0, 80,80, 0, 0, V]; 

%% Timing

%Days: simulation lenght
Days        = 400;
%timestep: for each day (e.g. 24 means every hour)
timestep    = 24;
%steps: simulation steps
ts          = 1;
steps       = 0:ts:Days*timestep-1; 

%Np: prediction horizen %days(?)
Np          = 5;       
%Nc: control horizen %days(?)
Nc          = 4;        

%% PBM model input

% Switching for paddlewheel on/off
swv         = 1;
% Reactor volume (L)
V           = 100;
% Influent characteristic % [O2, SS, SVFA, SIC, SH2, SIN, SIP, SI, XPB_ph, XPB_ch, XPB_an, XAHB, XAN, XS, XI];
Indata      = Influent';

% Changes in Indata
% Indata(3) = 1500;

% Light intensisty (W/m2)
intensity   = 54;                       

% Hydraulic retention time (day) 
HRT         = 4.0;
% Sludge retention time (day)
SRT         = 4.0;
% Area, m^2 (h = 0.20m -->20 cm )
A           = 0.5;
% Hight, m
h           = V/A/1000;
% HRT/SRT ratio defines the fraction of removed particles
fHS         = HRT/SRT;                  

% pw_b: Paddlewheel %If "1" paddlewheel is on
pw_b              = ones(24,1); 
pw_b(12:end)      = 1;                    

% timelight_b
timelight_b       = zeros(24,1);
% in case of natural lighting you may use
timelight_b(1:12) = 50/54*[15.6915261317682
                            26.6220146177290
                            41.0275523532562
                            57.4340255063066
                            73.0334577512867
                            84.3593891405176
                            88.5123655047883
                            84.3593891405176
                            73.0334577512867
                            57.4340255063066
                            41.0275523532562
                            26.6220146177290]; %intensity; 

% Metabolic switch
Ms_b              = ones(24,1)*0.28;
Ms_b(12:end)      = 0.28;

% Data organization
pw              = zeros(length(steps)+100,1);
light           = zeros(length(steps)+100,1);
Ms              = zeros(length(steps)+100,1);


for i = 0:(Days-1)

    Ms((1+(24*i)):(24*(i+1)))       = Ms_b;
    pw((1+(24*i)):(24*(i+1)))       = pw_b;
    mu = 55;
    sigma = 2;
    r = random('Normal',mu,sigma);
%     if i<40
        light((1+(24*i)):(24*(i+1)))    = [54*ones(12,1);zeros(12,1)]; 
%     else
%         light((1+(24*i)):(24*(i+1)))    = timelight_b;
%     end
end

% change in influent if needed
Indata(3) = 3000;
% assign an initial set-point for the controller
w         = 940*ones(length(steps),1);

% input/output flow rate
Qin  = zeros(length(steps),1);
Qout = zeros(length(steps),1);


%% ODE input and output variables

Input   = [light,pw,Ms];
Input_e = [];
options = odeset('NonNegative',1:16);
Output  = zeros(length(steps)+10,16);


%% Prediction model: ARMA, BJ, switched models, ...

%model order
% order of the system with respect to x (noise-free output)
na  = 1;
% order of the system with respect to u (input) 
nb  = 1;
% order of the system with respect to w (disturbance)
nc  = 1;
% order of the system with respect to v (input-free output)
nd  = 1;     
order = [na,nb,nc,nd];

% initial and reseting P 
P0 = 10e10*eye(na+nb+1);
P  = P0;

% initial guess for theta for adaptive identification
rng(120);
theta = rand(na+nb+1,1);          

%% Control Main Loop: simulate the plant and compute the control action
% defining an event when it happens, control action should be taken.
% if event = 1 --> system working in continuous mode.
% if event = 24 --> system working fed-batch based on a feeding schedule every 24h
% notice: it is assumed that ts = 1.
event = 24;
% starting point of the loop: can be assigned the user manually
ts1 = 1;        %event;     % na+nb+1;
k   = ts1 - 1;
% making u for further calculations
u     = 1*V/HRT*ones(length(steps),1);       % V/HRT*ones(floor(numel(steps)/event)+1,1);  
umax  = 2*V/HRT;
umin  = 0;

% numerating
s = 1;

% making yp and yc for further calculations
% yp     = 0*ones(length(steps),1);
% yc     = 0*ones(length(steps),1);
yd     = 0*ones(length(steps),1);
DeltaU_opt = 25*ones(Nc,1);

for i = ts1:event:numel(steps)    %numel(steps)-1

    % defining a helping variable to aviad index 0 for calling vector 'u'
    j = i + 1;

    % plant simulation
    % assign initial condition based on previous loop
    % This loop assign the first initial condition based the given initial.
    % the rest will be assigned based on outputs of the simulated plant.
    if k <= event    % k < na+nb+1
        k = k + event + 1;
    else
        initial = Output(i-1,:);                % initial = Output(i-1,:);
    end
    
    
    % tspan : time interval of simulation between step i and i+1
    % if we want to simulate the plant in a smaller timestep than 1hr i.e. ts ≠ 1, we can assign it here
    ts      = 1;
    tspan   = i-1:ts:i+event-2;       % steps(i):ts:steps(i+event);  % tspan   = [steps(i), steps(i+1)];        

    % assigning other inputs needed for the plant simulation
    % effective light, paddle wheel, and Ms
    
    light_e   = Input(i:i+event-1,1);                           % light_e   = Input(floor(i/ts)+1,1);
    pw_e      = Input(i:i+event-1,2);                           % pw_e      = Input(floor(i/ts)+1,1);
    Ms_e      = Input(i:i+event-1,3);                           % Ms_s      = Input(floor(i/ts)+1,1);

    % changes in operational conditions
    % influent VFA
    if i>2400
        Indata(3) = 3600;
    end
    
    % light conditions
    if i>4800 && i<=7200
        light_e = [60*ones(12,1);zeros(12,1)];
    elseif i>7200
        light_e = [50*ones(12,1);zeros(12,1)];
    end
    
    % metabolic switch constant
    if i>7200
        Ms_e = 0.32*ones(24,1);
    end


    Qin(i:i+event-1,1)  = [u(j-1);zeros(event-1,1)];            % Qin(floor(i/ts)+1,1)  = u(i-1);
    Qout(i:i+event-1,1) = [u(j-1);zeros(event-1,1)];            % Qout(floor(i/ts)+1,1) = u(i-1);


    Input_e   = [Input_e;light_e,Qin(i:i+event-1,1),Qout(i:i+event-1,1),pw_e,Ms_e]; 
    
    % simulate the plant between the tspan interval
   
    [~, yy]    = ode15s(@(t,y) PBM_MPC(t,y,Indata,Input_e,fHS,h), tspan, initial, options);
    % saving the output of the plant
    % [O2, SS, SVFA, SIC, SH2, SIN, SIP, SI, XPB_ph, XPB_ch, XPB_an, XAHB, XAN, XS, XI, V] 

    Output(i:i+event-1,:) = yy(1:1/ts:event/ts,:);
    %XPB = XPB_ph + XPB_ch + XPB_an --> output(:,9)+output(:,10)+output(:,11)
    yd(i:i+event-1) = sum(Output(i:i+event-1,[9 10 11]),2); 

    %%{

    %% adaptive system identification
    %
    % reseting P matrix to avoid approaching zero 
%     if (mod(i,50)== 0)     
%         P = P0;
%     end

    % making phi vector: the vector of the past inputs and outputs
    % to avoid empty phi

    if j-nb-1 <= 0
        if i+event-na*event <= 0
        phi = [-yp(i)*ones(na,1); u(j-1)*ones(nb+1,1)];
        else
            phi = [-yd(i:-event:i+event-na*event); u(j-1)*ones(nb+1,1)];
        end
    else
        phi = [-yd(i:-event:i+event-na*event); u(j-1:-1:j-(nb+1))];     % [-yp(i+event-1-1:-1:i+event-1-na); u(j-1:-1:j-(nb+1))];  or % phi = [-y(i-1:-1:i-na); u(i-1:-1:i-(nb+1))];
    end

    % call the ASI function
    [theta,P] = ASI(phi, yd(i+event-1,1), theta, P);
    % estimating yhat based on the identified system
    yhat(i)   = phi'*theta;

    % having the identified model i.e. A abd B matrices for MPC calculation
    A         = [1  theta(1:na)'];
    B         = [theta(na+1:end)']; 


    % call MPC analytical solution
    [G_u, Matrix_F, Matrix_G] = MPC_analytical_solution(A, B, Np , Nc);
    
    % calculate free response : f = G_u.∆.u + Matrix_F.Y
    % making ∆.u
    if i == 1
        DeltaU = zeros(nb,1);
    else
        DeltaU = u(j-1:-1:j-nb) - u(j-2:-1:j-nb-1);       % u(i-1:-1:i-(numel(B)-1))-u(i-2:-1:i-(numel(B)-1)-1);    % if nb>2, we may get error at j-(nb-1)-1 in first loop
    end
    % making Y
    if i+event-na*event-1 <= 0
        Y = [yd(i+event-1);yd(i)*ones(na,1)];
    else
        Y = yd(i+event-1:-event:i+event-na*event-1);       % y(i:-1:i-(numel(A)-1));
    end
    f      = G_u*DeltaU + Matrix_F*Y;

    % changing trajectory 1: quantity-driven | reduce outlet VFA
    if i>= 1*event
        if Output(i+event-1,3)>=250
            w(i+1:i+Np) = w(i+1-event:i+Np-event) + 5*ones(Np,1);
        elseif Output(i+event-1,3)<250 && Output(i+event-1,3)>=150
            w(i+1:i+Np) = w(i+1-event:i+Np-event);% + 5*ones(Np,1);
        elseif Output(i+event-1,3)<150
            w(i+1:i+Np) = w(i+1-event:i+Np-event) - 20*ones(Np,1);
        end
    end

    % changing trajectory2: quality-driven | increase outlet VFA
%     if i>= 5*event
%         if Output(i+event-1,3)/Indata(3)<=1/3
%             w(i+1:i+Np) = w(i+1-event:i+Np-event) - 10*ones(Np,1);
%         elseif Output(i+event-1,3)/Indata(3)>1/3 && Output(i+event-1,3)/Indata(3)<=2/5
%             w(i+1:i+Np) = w(i+1-event:i+Np-event); % - 5*ones(Np,1);
%         elseif Output(i+event-1,3)/Indata(3)>2/5
%             w(i+1:i+Np) = w(i+1-event:i+Np-event) + 10*ones(Np,1);
%         end
%     end
    ww(s) = w(i+1);
    s = s+1;

    % calculate control gain : K = (G'G + λI)^-1 G' 
    % then control actions within control prediction: DeltaU = K.(w - f) 
    % lambda as a controller design variable
    lambda = 100;
    % DeltaU = (Matrix_G'*Matrix_G + lambda*eye(size(Matrix_G,2)))\Matrix_G'*(w(i+1:i+Np) - f);    % ...*(w(i+1:i+Np) - f);
    % difference is only () for the last term
    DeltaU = (Matrix_G'*Matrix_G + lambda*eye(size(Matrix_G,2)))\(Matrix_G'*(w(i+1:i+Np) - f));    % ...*(w(i+1:i+Np) - f);
    %f(1)

%     Cost = @(u) (u'*(Matrix_G'*Matrix_G + lambda*eye(size(Matrix_G,2)))*u + 2*transpose(f - w(i+1:i+Np))*Matrix_G*u); % + transpose(f - w(i+1:i+Np))*(f - w(i+1:i+Np)));
% 
%     DeltaU_opt = fmincon(Cost,DeltaU_opt,[],[],[],[],0*ones(Nc,1),40*ones(Nc,1),[],[]);

    % so control action at time i based on DeltaU
    u(j:j+event-1) = u(j-1) + DeltaU(1); % DeltaU_opt(1);
    
%     u(j:j+event-1) = u(j-1) + DeltaU_opt(1);
    
    
    % saturation to avoid numerical problem

    if u(j) > 50
        u(j:j+event-1) = 50;
    elseif u(j) < umin
        u(j:j+event-1) = 0;
    else 
        u(j) = u(j);
    end

    %}

end

steps = 0:24:9600-1;
xb_control    = Output(:,9)+Output(:,10)+Output(:,11);

f1 = figure(13);
plot(steps, xb_control(1:24:9600,1),'k','linewidth',1.2)
hold on
plot(steps, ww,'r-s','MarkerSize',4)
ylabel('PPB (mg COD L$^{-1}$)','Interpreter','Latex','fontsize',12)
xlabel('Time (day)','Interpreter','Latex','fontsize',12)
xlim([0 9600])
ylim([700 1150])
yticks([700 915 980 1025 1110 1200])
yticklabels({'700', '915', '980', '1025', '1110', '1200'})
hold on 
legend ('Process response', 'Set-point','Outlet VFA','Interpreter','Latex','fontsize',10,'Location', 'east')
grid on
xticks([0 1200 2400 3600 4800 6000 7200 8400 9600])
xticklabels({'0', '50', '100', '150', '200', '250', '300', '350', '400'})
f1.Position = [100 100 1100 450];

f2 = figure(14);
plot(steps, Input_e(1:24:9600,2),'k','linewidth',1.2,'MarkerSize',5)
hold on 
plot(steps, 40*ones(400,1),'b-.','linewidth',1.2)
plot(steps, 0*ones(400,1),'b-.','linewidth',1.2)
legend ('Feeding flow rate','Upper bound', 'Lower bound', 'Interpreter','Latex','fontsize',10,'Location', 'northeast')
grid on
ylabel('u (L h$^{-1}$)','Interpreter','Latex','fontsize',12)
xlabel('Time (day)','Interpreter','Latex','fontsize',12)
xlim([0 9600])
ylim([0 40])
xticks([0 1200 2400 3600 4800 6000 7200 8400 9600])
xticklabels({'0', '50', '100', '150', '200', '250', '300', '350', '400'})
yticks([0 15 20 25 30 35 40])
f2.Position = [100 100 1100 450];
