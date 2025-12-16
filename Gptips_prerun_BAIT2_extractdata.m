%% this script is to extract data from BAIT2 for Millero equaions
clear;
ncopen('BAIT_2_1m_ptrc_Y2600_ETOPO_NPOFe.nc');% Nutrients, DO, Fer
ncopen('BAIT_2_1m_diad_Y2600_etopo_CFE_pH.nc'); % CFe
ncopen('NEMONEW_1m_T_S_ETOPO.nc'); % T, S

%% clear unreal value
PO4(PO4>1e10)=NaN;
NO3(NO3>1e10)=NaN;
O2(O2>1e10)=NaN;
SI(SI>1e10)=NaN;
FER(FER>1e10)=NaN;

CFE(CFE>1e10)=NaN;
PH(PH>1e10)=NaN;

SALIN(SALIN>1e10)=NaN;
TEMPER(TEMPER>1e10)=NaN;

%% select season
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

% deptht_gp=repelem(deptht, 360*180, 1);

% lon_gp = repelem(ETOPO60X, 180 * 31);  % 360*180*31 = 2,008,800×1
% lat_gp = repelem(ETOPO60Y, 360 * 31);  % 180*360*31 = 2,008,800×1

%potential
deptht_gp = repmat(reshape(deptht, 1, 1, []), 360, 180); % 3D matrix first
deptht_gp = deptht_gp(:); % then flatten

[lon2d, lat2d] = meshgrid(ETOPO60X, ETOPO60Y);  % 360x180
lon2d = lon2d';
lat2d = lat2d';
lon3d = repmat(lon2d, 1, 1, 31);
lat3d = repmat(lat2d, 1, 1, 31);

lon_gp = lon3d(:);
lat_gp = lat3d(:);
%% Fe solubility 
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
%% remove NaN of CFe
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9; % keep open ocean DFe and CFe/DFe between 10~90%

% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,14)>=60&BAIT2(:,14)<=90; % Arctic
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=20&BAIT2(:,13)<=120&BAIT2(:,14)>=-60&BAIT2(:,14)<=30; % Indian ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=280&BAIT2(:,13)<=380&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Atlantic ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,13)>=120&BAIT2(:,13)<=290&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Pacific ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=6000&BAIT2(:,14)<=-60; % Southern ocean

% % keep open ocean DFe and CFe/DFe between 10~90% and region

% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=3000&BAIT2(:,12)<=6000;
% % keep open ocean DFe and CFe/DFe between 10~90% and depth deep
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=200; %  surface
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=500&BAIT2(:,12)<=3000; % mesopelagic

% Upper ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500; % upper ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,14)>=60&BAIT2(:,14)<=90; % Arctic
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=20&BAIT2(:,13)<=120&BAIT2(:,14)>=-60&BAIT2(:,14)<=30; % Indian ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=280&BAIT2(:,13)<=380&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Atlantic ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,13)>=120&BAIT2(:,13)<=290&BAIT2(:,14)>=-60&BAIT2(:,14)<=60; % Pacific ocean
% rowskeep = BAIT2(:,8)>0&BAIT2(:,8)<=1&BAIT2(:,10)>=0.1&BAIT2(:,10)<=0.9&BAIT2(:,12)>=0&BAIT2(:,12)<=500&BAIT2(:,14)<=-60; % Southern ocean

BAIT2=BAIT2(rowskeep,:);
save("/Users/chengwangwang/Documents/MATLAB/CFemodel/Gptips_CFe_BAIT2_openocean_upper_Pacific.mat",'BAIT2')

%% stepwise regression
% C=B(:,[11:12,14,16:17,20,22]); % T, S, DO, DIP, Si, NO3, DFe
% C.zcfe=B.CFe_DFe_DOUBLE;
% mdl=stepwiselm(C,'ResponseVar','zcfe','Verbose',2)