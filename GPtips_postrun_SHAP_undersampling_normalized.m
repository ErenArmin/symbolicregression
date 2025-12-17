%% This script is to calculate SHAP of discovered equation from different undersampling scenarios.

clear
close all
%% Define symbolic variables and equations
syms tempis pH zsal;

% Define all 6 equations (replace with your actual equations)
% % Upper normalized
% eq1 = 0.00150822.*zsal - 0.08833078.*tempis - 0.03737602.*pH + 0.16679785;
% eq2 = 0.00987399.*tempis.^3 - 0.10449047.*tempis - 0.03755891.*pH + 0.16900027;
% eq3 = 0.00746506.*zsal - 0.12310007.*tempis - 0.03441699.*pH + 0.02563407.*exp(tempis) + 0.12747312;
% eq4 = 0.02015517.*exp(tempis) - 0.03557038.*pH - 0.12656072.*tempis + 0.00857834.*tempis.^3 + 0.13779068;
% eq5 = 0.01664604.*tempis.*pH - 0.02639142.*pH - 0.08519853.*tempis + 0.00864589.*tempis.^2 + 0.15239670;
% eq6 = 0.01555368.*tempis.^3 - 0.02640168.*pH - 0.00454439.*tempis.*exp(-1.0.*pH) - (1.0.*(9387298.0.*tempis.^2 + 4693649.0.*tempis + 9387298.0.*zsal))./(8388608.0.*tempis.*zsal.*pH - 1241002245.85196943) - 0.11371624.*tempis + 0.00037483.*pH.^3 + 0.15816562;
% 
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% % Deep normalized
% eq1 = 0.02499401.*pH.^2 - 0.12001890.*pH - 0.03395347.*tempis + 0.41679507;
% eq2 = -0.00372413.*tempis.^3 - 0.02265914.*tempis + 0.02566288.*pH.^2 - 0.12008120.*pH + 0.41917899;
% eq3 = 0.03320877.*exp(pH) - 0.14591706.*pH - 0.02479160.*tempis - 0.00309579.*tempis.^3 - 0.00904931.*pH.^3 + 0.39045256;
% eq4 = 0.03941372.*pH.^2 - 0.08231462.*pH - 0.02643130.*exp(pH) - 0.00281008.*tempis.^3 - 0.02683510.*tempis + 0.44581714;
% eq5 = 0.01115298.*tempis.^3 - 0.09648359.*pH - 0.02219198.*exp(tempis) - 0.01640956.*exp(pH) - 0.01640956.*tempis + 0.03397209.*pH.^2 + 0.46520537;
% eq6 = -0.01396910.*tempis.^2 - 0.02327966.*tempis + 0.02509216.*pH.^2 - 0.12132832.*pH + 0.00309626.*zsal + 0.43066709;
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% % Abyssal normalized
% eq1 = 0.68387759 - 0.08273422 .* pH;
% eq2 = 0.00515918 .* tempis - 0.08470383 .* pH + 0.68387735;
% eq3 = 0.00534886 .* pH.^2 - 0.08908364 .* pH + 0.67853129;
% eq4 = 0.00532860 .* pH.^2 - 0.08911766 .* pH + 0.00015253 .* tempis + 0.67855102;
% equations = {eq1, eq2, eq3, eq4};

% % Atlantic_terminate_normalized
% eq1 = 0.01755103.*tempis.^2 - 0.09291156.*tempis + 0.01835451.*pH.^2 - 0.06058839.*pH + 0.19687292;
% eq2 = 0.01757180.*tempis.^2 - 0.09584622.*tempis + 0.01815237.*pH.^2 - 0.06059993.*pH + 0.00322105.*zsal + 0.19705417;
% eq3 = 0.01069768.*tempis.^2 + 0.02139536.*tempis.*pH - 0.09650024.*tempis + 0.01067486.*pH.^2 - 0.05512233.*pH + 0.19786026;
% eq4 = 0.02166273.*tempis.*pH - 0.04600757.*pH - 0.00240802.*exp(2.0.*pH) - 0.09727695.*tempis + 0.01083136.*tempis.^2 + 0.01352884.*pH.^2 + 0.20107379;
% eq5 = 0.01395244.*tempis.^2 + 0.01514146.*tempis.*pH - 0.09615885.*tempis + 0.01291726.*pH.^2 - 0.05603993.*pH - 0.00012289.*zsal + 0.19632350;
% eq6 = 0.00819508.*tempis.*zsal - 0.06114791.*pH - 0.11099605.*tempis + 0.01278218.*tempis.^3 + 0.02131912.*pH.^2 + 0.20043548;

% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% % Pacific_terminate_normalized
% eq1 = 0.04807381.*tempis.*pH - 0.05084208.*pH - 0.10027424.*tempis - 0.01093820.*pH.^3 + 0.18944418;
% eq2 = 0.01994456.*pH - 0.10364011.*tempis - 0.06782110.*exp(pH) + 0.01118193.*(tempis + 2.0.*pH).^2 + 0.24396916;
% eq3 = 0.02204295.*exp(tempis) - 0.00210359.*zsal - 0.04664701.*pH - 0.13317253.*tempis + 0.04151107.*tempis.*pH - 0.01190331.*pH.^3 + 0.15641592;
% eq4 = 0.02459319.*exp(tempis) - 0.06393508.*pH - 0.13718683.*tempis + 0.01543784.*exp(pH) + 0.03741257.*tempis.*pH - 0.01092485.*pH.^3 + 0.13486117;
% eq5 = 0.00920113.*exp(2.0.*tempis) - 0.04305751.*pH - 0.00780957.*exp(tempis.^2) - 0.15373416.*tempis + 0.05381699.*tempis.*pH - 0.01093409.*pH.^3 + 0.17363650;
% eq6 = 1.08772365e-4.*exp(exp(tempis)) - 0.22888036.*tempis + 0.09575038.*exp(tempis) - 0.04699211.*exp(pH) + 0.01093028.*(tempis + 2.0.*pH).^2 - 0.02614765.*tempis.^4 + 0.10730731;
% 
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% % % Arctic_terminate_normalized
% eq1 = - 0.00218627 .* tempis.^3 + 0.00007020 .* zsal - 0.00752082 .* pH + 0.25317192;
% eq2 = - 0.00003711 .* pH.^3 - 0.00676292 .* pH + 0.00726522 .* tempis - 0.00186356 .* zsal - 0.00426395 .* exp(tempis) + 0.25854459;
% 
% equations = {eq1, eq2};
% 

% % % Indian_terminate_normalized
% eq1 = 0.03657123 .* pH.^2 - 0.06792670 .* pH - 0.08823091 .* tempis - 0.01440223 .* zsal + 0.20501813;
% eq2 = 0.03510053 .* pH.^2 - 0.06137312 .* pH - 0.15369749 .* tempis + 0.03670400 .* exp(tempis) + 0.14736268;
% eq3 = 0.05076465 .* exp(tempis) - 0.06168749 .* pH - 0.17147523 .* tempis - 0.01089541 .* tempis.^2 + 0.03538955 .* pH.^2 + 0.13531926;
% eq4 = 0.01378038 .* tempis.^2 + 0.02358734 .* tempis .* pH - 0.10478134 .* tempis + 0.02508776 .* pH.^2 - 0.05821639 .* pH - 0.00666808 .* zsal + 0.18810809;
% eq5 = 0.13049299 .* exp(tempis) - 0.06163472 .* pH - 0.26504263 .* tempis - 0.02785146 .* tempis.^4 + 0.03572610 .* pH.^2 + 0.04144833;
% eq6 = 0.00000003 .* tempis .* (1230723.0 .* pH - 5484824.0) - 0.05832560 .* pH + (0.00000012 .* zsal) ./ tempis + 0.02227873 .* tempis.^3 + 0.02245959 .* pH.^2 + 0.19518201;
% 
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% % % Southern_terminate_normalized
% eq1 = 0.43672967 - 0.10670976 .* pH;
% eq2 = 0.43672428 - 0.12446038 .* pH - 0.02091467 .* zsal;
% eq3 = 0.42905307 - 0.00904580 .* zsal .* pH - 0.11050622 .* pH;
% eq4 = 0.00331178 .* exp(pH) - 0.13153553 .* pH - 0.02576065 .* zsal + 0.00004630 .* (zsal - 1.0 .* pH).^3 + 0.01169612 .* tempis .* pH - 0.00331178 .* tempis.^2 + 0.42989942;
% eq5 = 0.01498577 .* pH.^2 - 0.00213357 .* exp(zsal) - 0.02588512 .* log(exp(2.0 .* pH)) - 0.05708447 .* pH + 0.42464893;
% 
% equations = {eq1, eq2, eq3, eq4, eq5};

% % IDP_terminate_normalized
% 
% eq1 = 0.03034265 .* zsal - 0.08885827 .* tempis - 0.09230398 .* pH + 0.26040432;
% eq2 = 0.02542491 .* exp(tempis) - 0.08257386 .* pH - 0.11353918 .* tempis + 0.21680341;
% eq3 = 0.04329578 .* exp(pH) - 0.11388072 .* pH - 0.08566749 .* tempis + 0.19741628;
% eq4 = 1.39506799 .* exp(tempis).^(1/2) - 0.08294655 .* pH - 0.80260465 .* tempis - 0.15661186 .* tempis.^2 - 0.00760032 .* tempis.^5 - 1.15924120;
% eq5 = 0.03034342 .* zsal - 0.16579218 .* tempis - 0.09048195 .* pH + 0.02233394 .* exp(tempis) + 0.03403940 .* exp(pH) + 0.17258251;
% eq6 = 0.01144050 .* exp(tempis) - 0.12108175 .* tempis - 1.0 .* pH .* (0.00241164 .* zsal - 0.03012688 .* tempis + 0.05176390) + (0.00001673 .* zsal) ./ pH + 0.00120582 .* zsal.^2 + 0.00120582 .* pH.^2 + 0.21589391;
% 
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

%% IDPDFe-IDP2025
eq1 = 0.025807.*pH.^2 - 0.089219.*pH - 0.082096.*tempis - 0.000722.*zsal + 0.251639;
eq2 = 0.011955.*tempis.^3 - 0.112365.*tempis + 0.027400.*pH.^2 - 0.082860.*pH + 0.242166;
eq3 = 0.013583.*tempis.^3 - 0.121425.*tempis + 0.027814.*pH.^2 - 0.081236.*pH + 0.007345.*zsal + 0.240678;
eq4 = 0.007300.*zsal - 0.110000.*tempis - 0.080000.*pH - 0.027000.*exp(tempis) + 0.030000.*tempis.^3 + 0.030000.*pH.^2 + 0.270000;
eq5 = -0.001933.*tempis.^6 + 0.035851.*tempis.^3 - 0.145653.*tempis + 0.028759.*pH.^2 - 0.080288.*pH + 0.005328.*zsal + 0.238867;
eq6 = 0.022780.*tempis.^3 - 0.023974.*tempis.^2 + 0.029090.*tempis.*pH - 0.135998.*tempis + 0.022970.*pH.^2 - 0.069182.*pH + 0.246580;

equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% IDPAtlantic
% eq1 = 0.26750000 - 0.09977000 .* pH - 0.07002000 .* tempis;
% eq2 = 0.02267293 .* tempis.^2 - 0.09767175 .* tempis - 0.08657492 .* pH + 0.24484059;
% eq3 = 0.01001063 .* (tempis + pH).^2 - 0.05622772 .* pH - 0.10734125 .* tempis + (0.00003427 .* zsal.^3 .* exp(zsal)) ./ pH + 0.23402822;
% eq4 = 0.00993621 .* zsal.^2 - 0.08106480 .* tempis - 0.09787738 .* pH + 0.25724474;
% eq5 = 0.00670519 .* zsal.^2 - 0.08544777 .* tempis - 0.12495459 .* pH + 0.03837655 .* exp(pH) + 0.20646805;
% eq6 = 0.00463372 .* pH .* (zsal + 2.0 .* pH) .* tempis.^2 - 0.09304525 .* tempis - 0.07248257 .* pH + 0.24997550;
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% IDPPacificIndian
% eq1 = 0.27013770 - 0.13167885 .* pH - 0.07749458 .* tempis;
% eq2 = 0.00433556 .* zsal - 0.08028272 .* tempis - 0.13223653 .* pH + 0.27013782;
% eq3 = 0.03559057 .* pH.^2 - 0.08464877 .* pH - 0.09620468 .* tempis + 0.23459780;
% eq4 = 0.06184726 .* tempis .* pH - 0.07609406 .* pH - 0.11834227 .* tempis + 0.22574866;
% eq5 = 0.03521567 .* pH.^2 - 0.07067921 .* pH - 0.16540030 .* tempis + 0.03643496 .* exp(tempis) + 0.17439157;
% eq6 = 0.03470098 .* pH.^2 - 0.07141653 .* pH - 0.17011568 .* tempis + 0.00417405 .* zsal + 0.03764347 .* exp(tempis) + 0.17289630;
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% IDPCruises - Pacific, Indian and Atlantic
% eq1 = 0.04286830 .* tempis .* pH - 0.07016125 .* pH - 0.12187321 .* tempis + 0.01469516 .* tempis.^2 + 0.22706600;
% eq2 = 0.01391207 .* (tempis + pH).^2 - 0.06423733 .* pH - 0.11710028 .* tempis + 0.22336136;
% eq3 = 0.02089748 .* (tempis + pH).^2 - 0.06535656 .* pH - 0.03053791 .* tempis .* pH - 0.11608139 .* tempis + 0.21976894;
% eq4 = 0.02086537 .* (tempis + pH).^2 - 0.00015294 .* zsal - 0.06532134 .* pH - 0.03037229 .* tempis .* pH - 0.11598653 .* tempis + 0.21976969;
% eq5 = 0.02212717 .* tempis.^2 + 0.01186122 .* tempis .* pH - 0.11730309 .* tempis + 0.02001843 .* pH.^2 - 0.06542579 .* pH + 0.21904069;
% eq6 = 0.01440451 .* tempis.^2 .* pH - 0.07375903 .* pH - 0.12739709 .* tempis + 0.01448409 .* (tempis + pH).^2 + 0.22020914;
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% IDPCruises - Southern ocean, Arctic, Pacific, Indian and Atlantic
% eq1 = 0.02584326 .* pH.^2 - 0.10427205 .* pH - 0.08140241 .* tempis + 0.28406978;
% eq2 = 0.02639133 .* pH.^2 - 0.10220953 .* pH - 0.08834302 .* tempis + 0.00751045 .* zsal + 0.28352156;
% eq3 = 0.00779397 .* tempis.^3 - 0.10436492 .* tempis + 0.02604294 .* pH.^2 - 0.10037967 .* pH + 0.27651015;
% eq4 = 0.07444499 .* exp(pH) - 0.19226724 .* pH - 0.07527430 .* exp(tempis) - 0.05346980 .* tempis + 0.05265947 .* tempis.^3 + 0.28946945;
% eq5 = 0.05621431 .* exp(pH) - 0.15501444 .* pH - 0.06673024 .* exp(tempis) - 0.06105242 .* tempis + 0.04756401 .* tempis.^3 - 0.00741590 .* pH.^3 + 0.29943326;
% eq6 = 0.01179840 .* zsal - 0.07961927 .* tempis - 0.09202656 .* pH - 0.06218696 .* exp(tempis) + 0.04630397 .* tempis.^3 + 0.03344095 .* pH.^2 + 0.34639257;

% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% Upper 500m
% eq1 = 0.02023373 .* tempis .* pH - 0.03644848 .* pH - 0.08583194 .* tempis + 0.17364119;
% eq2 = 0.02371891 .* exp(tempis) - 0.03809568 .* pH - 0.12080801 .* tempis + 0.01945134 .* zsal .* pH + 0.14035505;
% eq3 = 0.02087002 .* exp(tempis) - 0.01034426 .* zsal - 0.03617157 .* pH - 0.11235450 .* tempis + 0.02712945 .* zsal .* pH + 0.14362071;
% eq4 = 0.02636055 .* exp(tempis) - 0.02856348 .* pH - 0.12444948 .* tempis + 0.01821862 .* zsal .* pH + 0.00574531 .* pH.^2 + 0.13065983;
% eq5 = 0.01809761 .* exp(tempis) - 0.02663197 .* pH - 0.11005403 .* tempis + 0.01741853 .* tempis .* pH + 0.00578944 .* pH.^2 + 0.14043686;
% eq6 = 0.01099682 .* exp(tempis) - 0.02372146 .* pH - 0.12902358 .* tempis + 0.02140662 .* tempis .* pH + 0.01526868 .* tempis.^3 + 0.00597513 .* pH.^2 + 0.15108228;
% equations = {eq1, eq2, eq3, eq4, eq5, eq6};

% %% Upper Atlantic
% eq1 = 0.14748856 - 0.01923041 .* pH - 0.06941395 .* tempis;
% eq2 = -0.00382300 .* tempis .* zsal.^2 - 0.06361383 .* tempis - 0.01886254 .* pH + 0.14694583;
% eq3 = -0.00543666 .* zsal .* tempis.^2 - 0.06040931 .* tempis - 0.01857786 .* pH + 0.14614080;
% eq4 = 0.01785821 .* zsal - 0.10386992 .* tempis - 0.01846355 .* pH + 0.01785821 .* exp(tempis) - 0.00286356 .* exp(zsal) + 0.12450965;
% eq5 = 0.01308958 .* tempis.^2 - 0.06672472 .* tempis - 0.01936307 .* pH + 0.00321230 .* exp(pH) + 0.12987824;
% equations = {eq1, eq2, eq3, eq4, eq5};

% %% Upper Pacific
% eq1 = 0.03164355 .* tempis .* pH - 0.02947829 .* pH - 0.07214969 .* tempis + 0.14406008;
% eq2 = -0.01209000 .* tempis .* pH.^2 - 0.02572000 .* pH - 0.06551000 .* tempis + 0.14580000;
% eq3 = 0.02394258 .* exp(tempis) - 0.02666151 .* pH - 0.10521977 .* tempis + 0.03347331 .* tempis .* pH + 0.10561624;
% eq4 = -0.01217488 .* tempis .* pH.^2 - 0.02399451 .* pH - 0.09596051 .* tempis + 0.02222019 .* exp(tempis) + 0.11049345;
% eq5 = 0.00136181 .* pH - 0.10871946 .* tempis + 0.02650534 .* exp(tempis) - 0.01582711 .* exp(pH) + 0.02718569 .* tempis .* pH + 0.00953594 .* pH.^2 + 0.11467971;
% equations = {eq1, eq2, eq3, eq4, eq5};

% %% Uppper Indian
% eq1 = 0.14628893 - 0.01194913 .* zsal - 0.01646625 .* pH - 0.07576486 .* tempis;
% eq2 = 0.14526831 - 0.01016805 .* zsal - 0.01472219 .* pH - 0.00697574 .* tempis .* zsal .* pH - 0.07450174 .* tempis;
% eq3 = 0.03373002 .* exp(tempis) - 0.01633555 .* pH - 0.11877592 .* tempis + 0.09531239;
% eq4 = 0.01739676 .* tempis.^2 - 0.07364971 .* tempis - 0.00160251 .* zsal - 0.01573680 .* pH + 0.12887282;
% equations = {eq1, eq2, eq3, eq4};
% 
% %% Upper Atlantic
% eq1 = 0.14748856 - 0.01923041 .* pH - 0.06941395 .* tempis;
% eq2 = -0.00382300 .* tempis .* zsal.^2 - 0.06361383 .* tempis - 0.01886254 .* pH + 0.14694583;
% eq3 = -0.00543666 .* zsal .* tempis.^2 - 0.06040931 .* tempis - 0.01857786 .* pH + 0.14614080;
% eq4 = 0.01785821 .* zsal - 0.10386992 .* tempis - 0.01846355 .* pH + 0.01785821 .* exp(tempis) - 0.00286356 .* exp(zsal) + 0.12450965;
% eq5 = 0.01308958 .* tempis.^2 - 0.06672472 .* tempis - 0.01936307 .* pH + 0.00321230 .* exp(pH) + 0.12987824;
% equations = {eq1, eq2, eq3, eq5};

% %% Upper Arctic
% eq1 = 0.00019819 .* (tempis + pH + 3.09190000).^2 - 0.00649264 .* pH - 0.00226843 .* tempis.^3 + 0.25047213;
% equations = {eq1};

% %% Upper Southern
% eq1 = 0.37881228 - 0.05898453 .* pH;
% eq2 = 0.37881190 - 0.06149985 .* pH - 0.00298179 .* zsal;
% eq3 = -0.00011556 .* zsal.^3 + 0.00215291 .* zsal - 0.05823993 .* pH + 0.37846327;
% eq4 = -0.00002430 .* zsal.^3 + 0.00118298 .* pH.^3 - 0.06322958 .* pH + 0.37779266;
% eq5 = 0.00001687 .* tempis.^3 + 0.00013269 .* tempis + 0.00544130 .* pH.^2 - 0.06365097 .* pH - 0.00022987 .* zsal + 0.37333798;
% equations = {eq1,eq2,eq3,eq4,eq5};

% %% Upper IDP
% eq1 = 0.01409935 .* exp(tempis) - 0.06055856 .* tempis + 0.09582579;
% eq2 = 0.00605929 .* tempis.^2 - 0.03460713 .* tempis - 0.00912772 .* pH + 0.11122868;
% eq3 = 0.00743634 .* tempis .* zsal - 0.00880545 .* pH - 0.00021592 .* (zsal - 1.0 .* pH).^2 - 0.03743731 .* tempis + 0.10995593;
% equations = {eq1,eq2,eq3};

% %% Upper IDPAtlantic
% eq1 = 0.02651061 .* exp(tempis) - 0.09694221 .* tempis + 0.09098364;
% eq2 = -0.00414400 .* tempis.^3 + 0.01084000 .* tempis.^2 - 0.04111000 .* tempis - 0.01243000 .* zsal - 0.00414400 .* pH + 0.11990000;
% eq3 = -0.00438033 .* tempis.^2 .* pH + 0.01397631 .* tempis.^2 - 0.00056110 .* tempis .* pH.^2 - 0.05941850 .* tempis - 0.00438033 .* pH + 0.11747897;
% equations = {eq1,eq2,eq3};

% %% Upper IDPPacificIndian
% eq1 = 0.00613502 .* zsal.^2 + 0.01059113 .* pH .* zsal - 0.05736690 .* tempis - 0.02087942 .* pH + 0.11373258;
% eq2 = 0.01579374 .* tempis.^2 - 0.05386667 .* tempis - 0.00046050 .* zsal - 0.02051261 .* pH + 0.10719462;
% equations = {eq1,eq2};

% %% Upper IDPCruises
% eq1 = 0.01346438 .* tempis .* pH - 0.00962663 .* zsal - 0.01197650 .* pH - 0.06117750 .* tempis + 0.13050539;
% eq2 = 0.01562551 .* tempis .* zsal - 0.02193400 .* log(exp(pH)) - 0.06047089 .* tempis + 0.00000240 ./ zsal + 0.12294839;
% eq3 = 0.00058685 .* zsal.^3 + 0.00490839 .* zsal.^2 + 0.00131073 .* zsal - 0.10620517 .* tempis - 0.02089884 .* pH + 0.03233020 .* exp(tempis) + 0.07790195;
% eq4 = 0.01675545 .* tempis.^2 - 0.01376894 .* pH - 0.00170136 .* pH .* log(tempis.^2) - 0.06148391 .* tempis - 0.00144463 .* pH.^3 + 0.11443970;
% eq5 = 0.03263278 .* exp(tempis) - 0.00112311 .* zsal - 0.01211478 .* pH - 0.10229873 .* tempis + 0.00945986 .* zsal .* pH - 0.00184404 .* tempis.^2 - 0.00184404 .* tempis.^3 + 0.08126418;
% equations = {eq1,eq2,eq3,eq4,eq5};

% %% Upper IDPCruises-big
% eq1 = 0.16136700 - 0.00066298 .* zsal - 0.03280240 .* pH - 0.08339190 .* tempis;
% eq2 = 0.02729588 .* exp(tempis) - 0.03149299 .* pH - 0.11654128 .* tempis + 0.11950489;
% eq3 = 0.01047310 .* zsal - 0.13117787 .* tempis - 0.03092276 .* pH + 0.03260154 .* exp(tempis) + 0.11136810;
% eq4 = -0.00742328 .* pH .* zsal.^2 - 0.12608808 .* tempis - 0.02314490 .* pH + 0.03228040 .* exp(tempis) + 0.11466528;
% eq5 = 0.02633738 .* exp(tempis) - 0.02025077 .* pH - 0.11659816 .* tempis - 0.01525794 .* zsal.^2 .* pH + 0.00920959 .* zsal.^2 + 0.11753621;
% eq6 = 0.02641985 .* exp(tempis) - 0.00914705 .* zsal - 0.02063339 .* pH - 0.11352837 .* tempis + 0.02250772 .* zsal .* pH + 0.11603710;
% equations = {eq1,eq2,eq3,eq4,eq5,eq6};

% %% Upper IDP-DFe
% eq1 = 0.01458293 .* tempis.^3 - 0.11156377 .* tempis - 0.04747007 .* pH + 0.19021629;
% eq2 = -0.00208673 .* pH.^3 - 0.03649976 .* pH - 0.11226529 .* tempis + 0.01856386 .* exp(tempis) + 0.15897627;
% eq3 = 0.00065140 .* log(pH.^2) - 0.03639476 .* pH - 0.11217093 .* tempis + 0.01868894 .* exp(tempis) - 0.00203757 .* pH.^3 + 0.15966682;
% eq4 = 0.04961577 .* exp(pH).^(1 ./ 2) - 0.05827802 .* pH - 0.08522998 .* tempis + 0.01890477 .* tempis .* pH + 0.13059463;
% eq5 = 0.01114273 .* tempis.^3 .* pH - 0.08782591 .* tempis - 0.00238369 .* pH.^3 - 0.02563719 .* pH + 0.18242040;
% eq6 = 0.02066034 .* pH - 0.11430504 .* tempis - (15143141.0 .* pH) ./ (67108864.0 .* tempis + 33554432.0 .* pH + 331058092.44160000) + 0.01676658 .* tempis.^3 + 0.17822589;
% equations = {eq1,eq2,eq3,eq4,eq5,eq6};

n_equations = length(equations);

%% 2. Convert to Table-Compatible Functions
eq_functions = cell(1, n_equations);
for i = 1:n_equations
    % First create standard function
    func = matlabFunction(equations{i}, 'Vars', [tempis, zsal, pH]);
    
    % Wrap it to accept table input
    eq_functions{i} = @(tbl) func(tbl.tempis, tbl.zsal, tbl.pH);
end

%% 3. Load and Prepare Data (Optimized)
load("Gptips_CFe_BAIT2_openocean.mat"); % replace here by upper if testing upper ocean scenarios
% Replace normalize number
load('Gptips_prerun_BAIT2_fe3sol_openocean_IDP_DFe_IDP2025_normalized.mat'); % replace here with different means and sigmas from their scenarios.

% Create Latin Hypercube Sample for better coverage
n_samples = 1000;
rng(42); % For reproducibility

% Generate LHS sample in normalized space
lhs_sample = lhsdesign(n_samples, 3); % 3 variables
ranges = [min(BAIT2(:,[1,2,7])); max(BAIT2(:,[1,2,7]))]; % Physical ranges

% Scale LHS sample to physical ranges
X_phys = lhs_sample .* (ranges(2,:) - ranges(1,:)) + ranges(1,:);

% Normalize the sample
X_norm = (X_phys - mu_train) ./ sigma_train;

% Create input table
input_data = array2table(X_norm, 'VariableNames', {'tempis', 'zsal', 'pH'});
% load('shap_matlab_perturbate_input.mat')

%% 4. Compute SHAP Values (Optimized)
shap_results = cell(n_equations, 1);
shap_metrics = struct('mean_abs_shap', [], 'mean_shap', [], 'std_shap', [], 'mean_contrib',[]);
shap_contrib = cell(n_equations,1);
shap_values = cell(n_equations,1);
threshold = 1e-3; % Values below this will be set to zero
var_names = {'tempis', 'zsal', 'pH'}; % Must match your input_data variables

for eq_idx = 1:n_equations
    fprintf('Processing equation %d/%d...\n', eq_idx, n_equations);
    
    % Create explainer
        rng(42);
    explainer = shapley(eq_functions{eq_idx}, input_data);
    
    % Compute SHAP values in parallel
    shap_results{eq_idx} = fit(explainer, input_data);
    
    % Calculate metrics
    vals = shap_results{eq_idx}.Shapley.Value; %3 x1000

   % STRICT VARIABLE PRESENCE CHECK - NEW IMPROVED VERSION
    % Get the actual symbolic variables in the equation
    sym_vars = symvar(equations{eq_idx});
    present_vars = ismember(var_names, sym_vars);
    
    % Force SHAP to zero for absent variables
    vals(~present_vars, :) = 0;

      % vals(abs(vals) < threshold) = 0; % Set small values to zero
    shap_contrib{eq_idx} = 100 .* abs(vals')./sum(abs(vals'),2); %1000 X3
    shap_metrics(eq_idx).mean_contrib = mean(shap_contrib{eq_idx},1);
    shap_metrics(eq_idx).mean_abs_shap = mean(abs(vals'), 1);
    shap_metrics(eq_idx).mean_shap = mean(vals', 1);
    shap_metrics(eq_idx).std_shap = std(vals', 0, 1);
    fprintf('Mean absolute SHAP: %.4f, %.4f, %.4f\n', shap_metrics(eq_idx).mean_abs_shap);

    %save for EMD
    shap_values{eq_idx} = vals;
    mean_shap(eq_idx,:) = shap_metrics(eq_idx).mean_shap;
    var_shap(eq_idx,:) = shap_metrics(eq_idx).std_shap;
end

% Save for EMD
contributions = shap_contrib;

save('/Users/chengwangwang/Documents/MATLAB/SHAP_test/shap_matlab_values_IDP_DFe_IDP2025_normalized.mat','n_equations','shap_values','mean_shap','var_shap','contributions');
