%% This script is to extract variables of Millero equations from model output where have the same coordinates and depths from observation datasets
clear
%% BAIT2 input
ncopen('BAIT_2_1m_ptrc_Y2600_ETOPO_NPOFe.nc');% Nutrients, DO, Fer
ncopen('BAIT_2_1m_diad_Y2600_etopo_CFE_pH.nc'); % CFe
ncopen('NEMONEW_1m_T_S_ETOPO.nc'); % T, S

%% Remove unreal value
PO4(PO4>1e10)=NaN;
NO3(NO3>1e10)=NaN;
O2(O2>1e10)=NaN;
SI(SI>1e10)=NaN;
FER(FER>1e10)=NaN;

CFE(CFE>1e10)=NaN;
PH(PH>1e10)=NaN;

SALIN(SALIN>1e10)=NaN;
TEMPER(TEMPER>1e10)=NaN;

% Fe solubility 
fe3sol=Col_equation_PISCES_fe3sol(TEMPER,SALIN,PH);

%% Obs input
% Input DFe from GEOTRACES_IDP2025 
C = readtable('Gptips_DFe_IDP2025.csv');
C.CFe_nmol_kg_ = C.DFe_nmol_kg_ - C.SFe_nmol_kg_;
B = C;

%% transfer coords in observational dataset to match PISCES model - In PISCES, longitudes 20-380 (E->W), latitudes -90 (S) - 90 (N)   

for i=1:length(B.Longitude_degrees_east_)

    if B.Longitude_degrees_east_(i) < 20
        B.BAIT2_longitude(i) = 360 + B.Longitude_degrees_east_(i); % 1E ->361
    else B.BAIT2_longitude(i) = B.Longitude_degrees_east_(i);
    end
end

%% transfer coords and depth 
[~, B.BAIT2_longitude_index] = min(pdist2(B.BAIT2_longitude, ETOPO60X), [], 2);
[~, B.BAIT2_latitude_index] = min(pdist2(B.Latitude_degrees_north_, ETOPO60Y), [], 2);
[~, B.BAIT2_depth_idx] = min(pdist2(B.Depth_m_, deptht), [], 2);

%% extract data from BAIT2 referring coords and depths from observational datasets
for i=1:length(B.BAIT2_latitude_index)
    B.BAIT2_temperature (i) = TEMPER(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_salinity (i) = SALIN(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);

    B.BAIT2_DIP (i) = PO4(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_DIN (i) = NO3(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_Si (i) = SI(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_DO (i) = O2(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_PH (i) = PH(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);
    B.BAIT2_fe3sol (i) = fe3sol(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1);

    B.BAIT2_DFe (i) = FER(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1).*1e3;
    B.BAIT2_CFe (i) = CFE(B.BAIT2_longitude_index(i),B.BAIT2_latitude_index(i),B.BAIT2_depth_idx(i),1).*1e9;
    B.BAIT2_CFe_DFe = B.BAIT2_CFe ./ B.BAIT2_DFe;

    B.BAIT2_latitude (i) = ETOPO60Y(B.BAIT2_latitude_index(i));
    B.BAIT2_longitude (i) = ETOPO60X(B.BAIT2_longitude_index(i));
    B.BAIT2_depth (i) = deptht(B.BAIT2_depth_idx(i));

end

%% transfer longiotude back to real world (20-380 -> 0-360)
for i=1:length(B.Longitude_degrees_east_)

    if B.BAIT2_longitude_index(i) > 360
        B.BAIT2_longitude(i) = B.BAIT2_longitude(i) - 360; % 340 + y  = coords in PISCES_ETOPO if x shown on 360 scale
    end
end

%% Remove rows containing NaN
results = B;
results_cleaned = results(~isnan(results.BAIT2_CFe), :);

% data filteration
rowskeep = results_cleaned.BAIT2_DFe>0&results_cleaned.BAIT2_DFe<=1&...
results_cleaned.BAIT2_CFe_DFe>=0.1&results_cleaned.BAIT2_CFe_DFe<=0.9&...
results_cleaned.BAIT2_depth<=500;

% % data filteration for upper ocean
% rowskeep = results_cleaned.BAIT2_DFe>0&results_cleaned.BAIT2_DFe<=1&...
% results_cleaned.BAIT2_CFe_DFe>=0.1&results_cleaned.BAIT2_CFe_DFe<=0.9&...
% results_cleaned.BAIT2_depth<=500;

results_cleaned = results_cleaned(rowskeep,:);

%% save for GPTIPS
BAIT2_obs(:,1) = results_cleaned.BAIT2_temperature;
BAIT2_obs(:,2) = results_cleaned.BAIT2_salinity;
BAIT2_obs(:,3) = results_cleaned.BAIT2_DIP;
BAIT2_obs(:,4) = results_cleaned.BAIT2_DIN;
BAIT2_obs(:,5) = results_cleaned.BAIT2_Si;
BAIT2_obs(:,6) = results_cleaned.BAIT2_DO;
BAIT2_obs(:,7) = results_cleaned.BAIT2_PH;
BAIT2_obs(:,8) = results_cleaned.BAIT2_DFe;
BAIT2_obs(:,9) = results_cleaned.BAIT2_CFe;
BAIT2_obs(:,10) = results_cleaned.BAIT2_CFe_DFe;
BAIT2_obs(:,11) = results_cleaned.BAIT2_fe3sol;
BAIT2_obs(:,12) = results_cleaned.BAIT2_depth;
BAIT2_obs(:,13) = results_cleaned.BAIT2_longitude;
BAIT2_obs(:,14) = results_cleaned.BAIT2_latitude;

save("/Users/chengwangwang/Documents/MATLAB/CFemodel/Gptips_CFe_BAIT2_openocean_IDP_DFe_IDP2025.mat",'BAIT2_obs')
