
%"Str" defines the first hour over which simulation occurs
%"fin" defines the last hour of the simulation
str = 2905-1; 
fin = 6576; 

%Reading in VWC data to calculate R_soil bounds 
num = readmatrix('hourly2016');
VWC = num(str:fin,8);
rVWC = fillmissing(VWC,'linear');
for hn = 2:length(VWC)
    Rs(hn,1) = 0.1/(1.4018*rVWC(hn,1)+0.285);    
end

%BOUNDARIES OF SELECTED THERMOPHSYICAL PROPERTIES
%Note not all boundaries are defined here; some are defined in lines 36-37 
%and 42-43

%Soil thermal resistance
lbRs = min(Rs(2:length(Rs)));
ubRs = max(Rs(2:length(Rs)));

%Soil Thermal Capacitance in RC model
lbCsRC = 1101*724*0.1;
ubCsRC = 1101*873*0.1; 

%Volumetric heat capacity for canopy in FD model
lbCc = (790*1510); 
ubCc = (960*2248); 

%Volumetric heat capacity for soil in FD model
lbCsFD = 1101*724;
ubCsFD = 1101*873; 


%Estimate Boundaries for RC Model (Model 1):
%     Rcan   fveg   Csoil        Rsoil     LAI     Rins           Rconc          Cconc
lb = [0,     0,     lbCsRC,      lbRs,     0,      3.52*0.9,      0.2*0.8,       2080*800*0.3];     
ub = [3,     1,     ubCsRC,      ubRs,     3,      3.52*1.1,      0.2*1.2,       2400*1000*0.3];  


%Estimate Boundaies for FD Model (Model 2)
% %   Rcan   fveg   pCsoil       Rsoil     LAI    Rins          Rconc        pCconc
lb = [0,     0,     lbCsFD,      lbRs,     0,     3.52*0.9,     0.2*0.8,     2080*800];
ub = [3,     1,     ubCsFD,      ubRs,     3,     3.52*1.1,     0.2*1.2,     2400*1000];


%GENETIC ALGORITHM PROPERTIES
opts = optimoptions('ga','PlotFcn',@gaplotbestf,'CrossoverFrac',0.8,'PopulationSize',150,'StallGen',10,'Generations',10,'UseParallel',true);  
nvars = 8; % number of variables to optimize

%Only enable line 54 OR line 56 depending on the model being used
%RC Model (Model 1)
[x,fval,exitflag] = ga(@(x)RC_Model(x),nvars,[],[],[],[],lb,ub,[],[],opts);
 
 %FD Model (Model 2)
[y,fval,exitflag] = ga(@(y)FD_Model(y),nvars,[],[],[],[],lb,ub,[],[],opts);


