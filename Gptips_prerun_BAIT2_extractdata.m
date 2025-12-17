%% This script is to extract variables of Millero equations from the model outputs and subsampling them by different depths and coordinates. 
clear;
ncopen('BAIT_2_1m_ptrc_Y2600_ETOPO_NPOFe.nc');% Nutrients, DO, Fer
ncopen('BAIT_2_1m_diad_Y2600_etopo_CFE_pH.nc'); % CFe
ncopen('NEMONEW_1m_T_S_ETOPO.nc'); % T, S

%% Remove unreal value
PO4(PO4>1e10)=NaN;
NO3(NO3>1e10)=NaN;
O2(O2>1e10)=NaN;
SI(SI>1e10)=NaN;
FER(FER>1e10)=NaN; % dissolved Fe

CFE(CFE>1e10)=NaN; % colloidal Fe
PH(PH>1e10)=NaN;

SALIN(SALIN>1e10)=NaN;
TEMPER(TEMPER>1e10)=NaN;

%% Select season
i=1;
PO4_gp=PO4(:,:,:,i);
NO3_gp=NO3(:,:,:,i);
O2_gp=O2(:,:,:,i);
SI_gp=SI(:,:,:,i);
FER_gp=FER(:,:,:,i);

CFE_gp=CFE(:,:,:,i);
PH_gp=PH(:,:,:,i);

SALIN_gp=SALIN(:,:,:,i);
TEMPER_gp=TEMPER(:,:,:,i);

% input depth with every varaible
deptht_gp = repmat(reshape(deptht, 1, 1, []), 360, 180); % 3D matrix first
deptht_gp = deptht_gp(:); % then flatten

% input longitude with latitude (360*180)
[lon2d, lat2d] = meshgrid(ETOPO60X, ETOPO60Y);  % 360x180
lon2d = lon2d';
lat2d = lat2d';
lon3d = repmat(lon2d, 1, 1, 31);
lat3d = repmat(lat2d, 1, 1, 31);

lon_gp = lon3d(:);
lat_gp = lat3d(:);

% input Fe solubility 
fe3sol=Col_equation_PISCES_fe3sol(TEMPER_gp,SALIN_gp,PH_gp);


%% Compile
BAIT2(:,1) = TEMPER_gp(:);
BAIT2(:,2) = SALIN_gp(:);
BAIT2(:,3) = PO4_gp(:);
BAIT2(:,4) = NO3_gp(:);
BAIT2(:,5) = SI_gp(:);
BAIT2(:,6) = O2_gp(:);
BAIT2(:,7) = PH_gp(:);
BAIT2(:,8) = FER_gp(:).*1e3;
BAIT2(:,9) = CFE_gp(:).*1e9;
BAIT2(:,10) = (CFE_gp(:).*1e9)./(FER_gp(:).*1e3);
BAIT2(:,11) = fe3sol(:);
BAIT2(:,12) = deptht_gp;
BAIT2(:,13) = lon_gp;
BAIT2(:,14) = lat_gp;
%% Filtered by DFe 0 - 1 nmol/L to focus on the open ocean Fe cycling away from strong point sources
% each 'rowskeep' represent each scenario (e.g., data only from Pacific or data only from surface)

rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9; % Ideal scenario with global data 
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,14)>=60&BAIT2(:,14)<=90; % Arctic
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=20&BAIT2(:,13)<=120&BAIT2(:,14)>=-60&BAIT2(:,14)<=30; % Indian ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=280&BAIT2(:,13)<=380&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Atlantic ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=120&BAIT2(:,13)<=290&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Pacific ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,14)<=-60; % Southern ocean

% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=200; % upper 
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=500&BAIT2(:,12)<=3000;% deep
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=3000&BAIT2(:,12)<=6000; % Abyssal

% Upper ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500; % upper ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,14)>=60&BAIT2(:,14)<=90; % upper ocean in Arctic
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=20&BAIT2(:,13)<=120&BAIT2(:,14)>=-60&BAIT2(:,14)<=30; % upper ocean in Indian ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=280&BAIT2(:,13)<=380&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % upper ocean at Atlantic ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=120&BAIT2(:,13)<=290&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % upper ocean at Pacific ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,14)<=-60; % Southern ocean

BAIT2=BAIT2(rowskeep,:);
save("/Users/chengwangwang/Documents/MATLAB/CFemodel/Gptips_CFe_BAIT2_openocean.mat",'BAIT2')
