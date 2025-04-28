%% 1. Land Earth System Model (LESTM) V1.0
%==========================================================================
% Lead Researchers:
%   Researcher Minming Cui
%   Researcher Xinbin Feng
%   Researcher Ping Li
%
% Usage:
%   This code is provided for academic research and teaching purposes only.
%   For any other use, or for redistribution and further development,
%   please contact Researcher Minming Cui at the Institute of Geochemistry,
%   Chinese Academy of Sciences, to obtain formal authorization.
%==========================================================================


clear;
start_year =;
end_year =;
years = start_year:end_year;
numYears = length(years);

demFile = 'your_data';
[DEM, R_dem] = readgeoraster(demFile, 'OutputType', 'double');
maskofElevation = DEM > 0;
DEM_initial = DEM;
soilThicknessFile = 'your_data';
soilFile = 'your_data';
numYears = length(years);

nodataValue = -9999;
[soil_thickness, grid, slope, aspect] = ...
    initSpatialParams(DEM, R_dem, soilThicknessFile, nodataValue, numYears);

Precipitation_mm_year = NaN(numRows, numCols, numYears);  
Evaporation_mm_year = NaN(numRows, numCols, numYears);  
Baseflow_kg_m2_year = NaN(numRows, numCols, numYears);    

scale_factor = 1;

for y = 1:numYears
    year = years(y);
    precipFile = sprintf('your_data', year);
    [precip_data, ~] = readgeoraster(precipFile, 'OutputType', 'double');
    precip_data(precip_data == -9999) = NaN;  
    Precipitation_mm_year(:,:,y) = precip_data * scale_factor * 365; 
    Precipitation_mm_year(:,:,y) = fillmissing(Precipitation_mm_year(:,:,y), 'constant', 0);
end


for y = 1:numYears
    year = years(y);
    evapFile = sprintf('your_data', year);
    [evap_data, ~] = readgeoraster(evapFile, 'OutputType', 'double');
    evap_data(evap_data == -9999) = NaN; 
    Evaporation_mm_year(:,:,y) = evap_data * scale_factor; 
    Evaporation_mm_year(:,:,y) = fillmissing(Evaporation_mm_year(:,:,y), 'constant', 0);
end

for y = 1:numYears
    year = years(y);
    baseflowFile = sprintf('your_data', year);
    [baseflow_data, ~] = readgeoraster(baseflowFile, 'OutputType', 'double');
    baseflow_data(baseflow_data == -9999) = NaN;  
    Baseflow_kg_m2_year(:,:,y) = baseflow_data * scale_factor; 
    Baseflow_kg_m2_year(:,:,y) = fillmissing(Baseflow_kg_m2_year(:,:,y), 'constant', 0);
end

if any(isnan(Precipitation_mm_year(:))) || any(isnan(Evaporation_mm_year(:))) || any(isnan(Baseflow_kg_m2_year(:)))
    Precipitation_mm_year(isnan(Precipitation_mm_year)) = 0;
    Evaporation_mm_year(isnan(Evaporation_mm_year)) = 0;
    Baseflow_kg_m2_year(isnan(Baseflow_kg_m2_year)) = 0;
end


soil_density = ; 
water_density = ; 

rock_density = ; 
k_bedrock = ;     
k_transport = ; 
m_sp = ;  
n_sp = ;  
sp_crit = ;
F_f = ;
solver = 'basic';
dt_min = ;  

[~, linearIdx] = min(DEM(:),[],'omitnan');
[xLow, yLow] = ind2sub(size(DEM), linearIdx);
drain = false(numRows, numCols);
drain(xLow, yLow) = true;
pool_matrix = false(numRows, numCols);
drainage = false(numRows, numCols);
min_elevation = 1e-3; 
validMask = DEM > 0;
S_initial = soil_thickness;  
S_initial(~maskofElevation) = NaN;  
S_initial(isnan(DEM)) = NaN;       
S = S_initial; 
Hbedrock_initial = DEM_initial - S_initial;  
Hbedrock_initial(~maskofElevation) = NaN;  
Hbedrock = Hbedrock_initial; 
hstar = ;   
A0 = ;      
A1 = ;    
b = ;   

useChemicalWeathering = true;
phi = ;         
Rn_max = ;     
f_w = ;         
C_eq = ;      
alpha = ;       
USP_initial_depth = ;  
USP_initial = USP_initial_depth * soil_density * pixel_area_ha .* maskofElevation;  
USP_initial(~isfinite(S_initial)) = NaN;  
USP = USP_initial; 

M_C = 12.01;
M_N = 14.01;
M_P = 30.97;
M_S = 32.07;
M_O = 16.00;

SOC_initial_value_t_ha = ;  
SIC_initial_value_t_ha = ; 
DON_initial_value_t_ha = ;  
DIN_initial_value_t_ha = ;  
POP_initial_value_t_ha = ; 
Pi_initial_value_t_ha = ; 
Sorg_initial_value_t_ha = ; 
Sinorg_initial_value_t_ha = ;  
O2_initial = ;    
OR_initial_t_ha = ;   
Cd_initial_mg_kg = ;  
Hg_initial_mg_kg = ;  
As_initial_mg_kg = ;   

soil_mass_t_ha_matrix = S .* pixel_area_ha .* soil_density;  
soil_mass_t_ha_matrix(isnan(soil_mass_t_ha_matrix) | soil_mass_t_ha_matrix == 0) = 1e-10;  
Cd_initial_mg_ha = Cd_initial_mg_kg .* soil_mass_t_ha_matrix * 1000; 
Hg_initial_mg_ha = Hg_initial_mg_kg .* soil_mass_t_ha_matrix * 1000;  
As_initial_mg_ha = As_initial_mg_kg .* soil_mass_t_ha_matrix * 1000;  
SOC = SOC_initial_value_t_ha * exp(-0.02 * slope) .* maskofElevation;   
SIC = SIC_initial_value_t_ha * exp(-0.01 * slope) .* maskofElevation;    
DON = DON_initial_value_t_ha * exp(-0.03 * slope) .* maskofElevation;   
DIN = DIN_initial_value_t_ha * exp(-0.02 * slope) .* maskofElevation;    
POP = POP_initial_value_t_ha * exp(-0.04 * slope) .* maskofElevation;   
Pi = Pi_initial_value_t_ha * exp(-0.03 * slope) .* maskofElevation;     
Sorg = Sorg_initial_value_t_ha * exp(-0.025 * slope) .* maskofElevation;  
Sinorg = Sinorg_initial_value_t_ha * exp(-0.02 * slope) .* maskofElevation;
O2 = O2_initial * ones(numRows, numCols) .* maskofElevation;     
OR = OR_initial_t_ha * ones(numRows, numCols) .* maskofElevation; 
Cd = Cd_initial_mg_ha .* maskofElevation;  
Hg = Hg_initial_mg_ha .* maskofElevation;  
As = As_initial_mg_ha .* maskofElevation;  

metal_generated_ratio_Cd = ;  
metal_generated_ratio_Hg = ; 
metal_generated_ratio_As = ;  
metal_eroded_ratio_Cd = ;  
metal_eroded_ratio_Hg = ;  
metal_eroded_ratio_As = ;  


soil_generated_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3); 
soil_eroded_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);    
soil_eroded_water_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);  
soil_eroded_wind_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);  

SOC_history_ha_year = zeros(numYears, numRows, numCols);  
SIC_history_ha_year = zeros(numYears, numRows, numCols);  

DON_history_ha_year = zeros(numYears, numRows, numCols);  
DIN_history_ha_year = zeros(numYears, numRows, numCols);  

POP_history_ha_year = zeros(numYears, numRows, numCols);  
Pi_history_ha_year = zeros(numYears, numRows, numCols);   

Sorg_history_ha_year = zeros(numYears, numRows, numCols);     
Sinorg_history_ha_year = zeros(numYears, numRows, numCols);  

O2_history_mg_L_year = zeros(numYears, numRows, numCols);      
OR_history_mol_O_ha_year = zeros(numYears, numRows, numCols);  

Cd_history_mg_ha_year = zeros(numYears, numRows, numCols);  
Hg_history_mg_ha_year = zeros(numYears, numRows, numCols); 
As_history_mg_ha_year = zeros(numYears, numRows, numCols);  

% [Cd, Hg, As]
soil_generated_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3); 
soil_eroded_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);   
soil_eroded_water_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);  
soil_eroded_wind_metals_history_mg_ha_year = zeros(numYears, numRows, numCols, 3);   
soil_generated_history_ha_year = zeros(numYears, numRows, numCols);  
soil_eroded_history_ha_year = zeros(numYears, numRows, numCols);    
soil_generated_total_history_ton_year = zeros(numYears, 1); 
soil_eroded_total_history_ton_year = zeros(numYears, 1);     


CNPS_water_history_ton_year = zeros(numYears, numRows, numCols);    
Metals_water_history_mg_year = zeros(numYears, numRows, numCols);   
E_water_history_ha_year = zeros(numYears, numRows, numCols);        
E_wind_history_ha_year = zeros(numYears, numRows, numCols);         
E_total_history_ha_year = zeros(numYears, numRows, numCols);        
D_history_ha_year = zeros(numYears, numRows, numCols);             

total_generated_soil_ton_year = ;      
total_water_eroded_ton_year = ;       
total_wind_eroded_ton_year = ;        
total_soil_eroded_ton_year = ;        
total_soil_eroded_underground_ton_year = ;  
total_water_inflow_ton_year = ;       
total_water_outflow_ton_year = ;       
soil_generated_total_history_ton_year = zeros(numYears, 1);  
soil_eroded_underground_history_ha_year = zeros(numYears, numRows, numCols); 

DIN_sim_ton_year = zeros(numYears, 1);   
Pi_sim_ton_year = zeros(numYears, 1);    
Sinorg_sim_ton_year = zeros(numYears, 1); 

Cd_sim_mg_year = zeros(numYears, 1);    
Hg_sim_mg_year = zeros(numYears, 1);    
As_sim_mg_year = zeros(numYears, 1);    

oceanLevelParameter = ;  

[ ocean, VOcean_m3, AOcean_m2, ZBeachLevel, border_mask, laplacianKernel ] = ...
    detectOceanAndKernel( DEM, pixel_area, oceanLevelParameter );


numVerticalLayers = ;  
layerThickness = ;   
total_depth = layerThickness * numVerticalLayers;  

diffusion_Cd = ;      
diffusion_Hg = ;      
diffusion_As = ;      
delta_z      = layerThickness;
dt           = ;

Cd_concentration_3D = zeros(numYears, numVerticalLayers, numRows, numCols);
Hg_concentration_3D = zeros(numYears, numVerticalLayers, numRows, numCols);
As_concentration_3D = zeros(numYears, numVerticalLayers, numRows, numCols);

Cd_initial = Cd;          
Hg_initial = Hg;         
As_initial = As;          

metalParams = struct();
metalParams.Vmax = struct('Cd', 200, 'Hg', 150, 'As', 180);          
metalParams.Km   = struct('Cd',  50, 'Hg',  30, 'As',  40);           
metalParams.leach_rate    = struct('Cd', 0.02, 'Hg', 0.01, 'As', 0.03);
metalParams.erosion_ratio = struct('Cd', 0.05, 'Hg', 0.03, 'As', 0.07);

cnpsParams = struct();
cnpsParams.generalRate = ;    
cnpsParams.epsilon     = ;   

max_Soi_mm =;      
max_Epi_mm = ;      
k_inf = ;          
k_rech = ;          
k_AET = ;          

Soi_mm = max_Soi_mm / 2 * ones(numRows, numCols) .* maskofElevation;   
Epi_mm = max_Epi_mm / 2 * ones(numRows, numCols) .* maskofElevation;   

Soi_history_mm_year = zeros(numYears, numRows, numCols);   
Epi_history_mm_year = zeros(numYears, numRows, numCols);   
Inf_history_mm_year = zeros(numYears, numRows, numCols);  

Eact_history_mm_year = zeros(numYears, numRows, numCols);  
Exc_Soi_history_mm_year = zeros(numYears, numRows, numCols);  
Exc_Epi_history_mm_year = zeros(numYears, numRows, numCols); 

Rech_history_mm_year = zeros(numYears, numRows, numCols);  
InSurf_history_mm_year = zeros(numYears, numRows, numCols); 
RunOff_history_mm_year = zeros(numYears, numRows, numCols); 

for t_step = 1:numYears
    current_precipitation_mm_year = Precipitation_mm_year(:,:,t_step);  
    current_pet_mm_year = Evaporation_mm_year(:,:,t_step);               
    current_baseflow_kg_m2_year = Baseflow_kg_m2_year(:,:,t_step);     

    precipitation_volume_L_year = current_precipitation_mm_year .* pixel_area / 1000 .* 1e3; 
    Cd_precip_input_mg_year = 0.0001 * precipitation_volume_L_year; 
    Hg_precip_input_mg_year = 0.00005 * precipitation_volume_L_year; 
    As_precip_input_mg_year = 0.0002 * precipitation_volume_L_year; 
    baseflow_volume_ton_year = (current_baseflow_kg_m2_year .* pixel_area) / (water_density * 1e3);  
    
    Cd_baseflow_input_mg_year = 0.001 * baseflow_volume_ton_year;  
    Hg_baseflow_input_mg_year = 0.0005 * baseflow_volume_ton_year; 
    As_baseflow_input_mg_year = 0.002 * baseflow_volume_ton_year;  

    Runoff_mm_year = current_precipitation_mm_year - current_pet_mm_year; 
    Runoff_mm_year(Runoff_mm_year < 0) = 0;  
    Runoff_mm_year(Runoff_mm_year > 1500) = 1500;
    Runoff_volume_m3_year = Runoff_mm_year .* pixel_area / 1000;  
    Runoff_mass_ton_year = Runoff_volume_m3_year * water_density;    

    Q = Runoff_volume_m3_year;  
    E_water_ha_year = k_bedrock * (Q).^m_sp .* (slope).^n_sp;  
    E_water_ha_year = E_water_ha_year - sp_crit;  
    E_water_ha_year(E_water_ha_year < 0) = 0;

    sand_fraction = ;          
    clay_fraction = ;       
    organic_matter = ;       
    surface_crust = ;          
    surface_roughness = ;        
    z0 = ;                     
    u_t = ;                   
    U_2 = ;                      
    rho = ;                   
    g = ;                    
    N = ;                        
    
    WF = (rho / g) * (1 / N) * sum((U_2 - u_t).^2 * U_2);

    EF = 1 - surface_crust;

    SCF = 1 / (1 + 0.0066 * clay_fraction + 0.021 * organic_matter)^2;
    K_prime = surface_roughness / z0;
    SF = sin(slope * (pi / 180)); 

    Q_t_max_RWEEQ = 109.8 * WF * EF * SCF * K_prime * SF;
    E_wind_ha_year = Q_t_max_RWEEQ .* maskofElevation; 
    H_initial = S;  

    f_solute = ;  
    E=E_water_ha_year;
    H_mc = log((E/1000 * rock_density)/10^4.01) / -2300 * 1000;

    q_local = Runoff_volume_m3_year ./ pixel_area;   
    Dw = calculate_Dw(H_mc, E_water_ha_year, phi, Rn_max, f_w, C_eq);
    Q_flux = calculate_Q(Dw, q_local, alpha);

    physical_wx = min(A1 .* exp(-H_initial / hstar), A0 + b .* H_initial);  
    total_weathering = physical_wx + Q_flux; 
    soil_wx = total_weathering * (1 - f_solute);    
    Water_SiO2_flux = total_weathering * f_solute;   

    soil_generated_ha_year = (rock_density / soil_density) * physical_wx * pixel_area_ha;  
    soil_generated_Cd = soil_generated_ha_year .* metal_generated_ratio_Cd * 1000;  
    soil_generated_Hg = soil_generated_ha_year .* metal_generated_ratio_Hg * 1000;  
    soil_generated_As = soil_generated_ha_year .* metal_generated_ratio_As * 1000;  

    soil_generated_metals_history_mg_ha_year(t_step, :, :, 1) = soil_generated_Cd;  
    soil_generated_metals_history_mg_ha_year(t_step, :, :, 2) = soil_generated_Hg; 
    soil_generated_metals_history_mg_ha_year(t_step, :, :, 3) = soil_generated_As; 

    Hbedrock = Hbedrock - total_weathering;  
    S = S + (rock_density / soil_density) * soil_wx;  
    Hbedrock(Hbedrock < min_elevation) = min_elevation;
    D_ha_year = (k_bedrock / k_transport) * E_water_ha_year ./ (Q + 1e-10);  

    total_erosion_ha_year = E_water_ha_year + E_wind_ha_year;  
    soil_eroded_Cd = total_erosion_ha_year .* metal_eroded_ratio_Cd * 1000;  
    soil_eroded_Hg = total_erosion_ha_year .* metal_eroded_ratio_Hg * 1000;  
    soil_eroded_As = total_erosion_ha_year .* metal_eroded_ratio_As * 1000;  

    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 1) = soil_eroded_Cd;  
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 2) = soil_eroded_Hg; 
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 3) = soil_eroded_As;  
    soil_eroded_water_mg_ha_year = E_water_ha_year .* metal_eroded_ratio_Cd * 1000;  
    soil_eroded_wind_mg_ha_year = E_wind_ha_year .* metal_eroded_ratio_Cd * 1000;    

    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 1) = E_water_ha_year .* metal_eroded_ratio_Cd * 1000;
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 2) = E_water_ha_year .* metal_eroded_ratio_Hg * 1000;  
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 3) = E_water_ha_year .* metal_eroded_ratio_As * 1000; 
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 1) = E_wind_ha_year .* metal_eroded_ratio_Cd * 1000; 
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 2) = E_wind_ha_year .* metal_eroded_ratio_Hg * 1000; 
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 3) = E_wind_ha_year .* metal_eroded_ratio_As * 1000;  

    delta_S_erosion_m = total_erosion_ha_year / (soil_density * 1e4);  
    delta_S_deposition_m = D_ha_year / (soil_density * 1e4);  

    S = S - delta_S_erosion_m + delta_S_deposition_m;  

    DEM = Hbedrock + S;  

    soil_generated_history_ha_year(t_step, :, :) = soil_generated_ha_year;  
    soil_eroded_history_ha_year(t_step, :, :) = total_erosion_ha_year;     
    soil_generated_total = nansum(soil_generated_ha_year(:) .* pixel_area_ha);  
    soil_eroded_total = nansum(total_erosion_ha_year(:) .* pixel_area_ha);      
    soil_generated_total_history_ton_year(t_step) = soil_generated_total;
    soil_eroded_total_history_ton_year(t_step) = soil_eroded_total;
    soil_net_change_ha_year = soil_generated_ha_year - total_erosion_ha_year;  

    k_percolation = 0.05; 
    soil_percolation_ton_year = k_percolation * S .* soil_density .* pixel_area_ha;  
    max_percolation_ton_year = S .* soil_density .* pixel_area_ha; 
    soil_percolation_ton_year = min(soil_percolation_ton_year, max_percolation_ton_year); 
    soil_percolation_m_year = soil_percolation_ton_year / (soil_density * pixel_area_ha);  

    S = S - soil_percolation_m_year;
    USP = USP + soil_percolation_ton_year;  
    k_underground_erosion = 0.02;  
    soil_eroded_underground_ha_year = k_underground_erosion * USP;  
    soil_eroded_underground_ha_year(~maskofElevation) = 0;  
    total_soil_eroded_underground_ton_year = total_soil_eroded_underground_ton_year + sum(soil_eroded_underground_ha_year(:) .* pixel_area_ha / (numRows * numCols)); 
    soil_eroded_underground_history_ha_year(t_step, :, :) = soil_eroded_underground_ha_year;  

    soil_eroded_underground_Cd = soil_eroded_underground_ha_year .* metal_eroded_ratio_Cd * 1000;
    soil_eroded_underground_Hg = soil_eroded_underground_ha_year .* metal_eroded_ratio_Hg * 1000; 
    soil_eroded_underground_As = soil_eroded_underground_ha_year .* metal_eroded_ratio_As * 1000; 


    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 1) = soil_eroded_underground_Cd;
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 2) = soil_eroded_underground_Hg;  
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 3) = soil_eroded_underground_As; 

    USP = USP - soil_eroded_underground_ha_year;  
    underground_loss_ha_year = soil_percolation_ton_year + soil_eroded_underground_ha_year; 
    underground_loss_history_ton_year(t_step, :, :) = underground_loss_ha_year;

    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 1) = E_water_ha_year .* metal_eroded_ratio_Cd * 1000;  
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 2) = E_water_ha_year .* metal_eroded_ratio_Hg * 1000;  
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 3) = E_water_ha_year .* metal_eroded_ratio_As * 1000;  

    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 1) = E_wind_ha_year .* metal_eroded_ratio_Cd * 1000;   
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 2) = E_wind_ha_year .* metal_eroded_ratio_Hg * 1000;   
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 3) = E_wind_ha_year .* metal_eroded_ratio_As * 1000;   

    Inf = k_inf * current_precipitation_mm_year; 
    Eact = k_AET * current_pet_mm_year; 
    Soi_new = Soi_mm + Inf - Eact;  

    if Soi_new > max_Soi_mm
        Soi_mm = max_Soi_mm;
        Exc_Soi = Soi_new - max_Soi_mm;
    elseif Soi_new >= 0
        Soi_mm = Soi_new;
        Exc_Soi = 0;
    else
        Soi_mm = 0;
        Exc_Soi = 0;  
        InSurf = 0 - Soi_new;  
    end

    Epi_new = Epi_mm + squeeze(Rech_history_mm_year(max(t_step-1,1), :, :)) .* 0.1 - 0.05 * Epi_mm; 
    Epi_new = max(Epi_new, 0); 
    Exc_Epi = Epi_new - max_Epi_mm;  
    Exc_Epi(Exc_Epi < 0) = 0;  
    Epi_mm = min(Epi_new, max_Epi_mm); 

    Rech = k_rech * Exc_Soi;  
    InSurf = Exc_Soi - Rech;  
    InSurf(InSurf < 0) = 0;  
    Rech_history_mm_year(t_step, :, :) = Rech;  

    Soi_history_mm_year(t_step, :, :) = Soi_mm;    
    Epi_history_mm_year(t_step, :, :) = Epi_mm;    
    Inf_history_mm_year(t_step, :, :) = Inf;        
    Eact_history_mm_year(t_step, :, :) = Eact;      
    Exc_Soi_history_mm_year(t_step, :, :) = Exc_Soi; 
    Exc_Epi_history_mm_year(t_step, :, :) = Exc_Epi;  
    InSurf_history_mm_year(t_step, :, :) = InSurf;    
    RunOff_history_mm_year(t_step, :, :) = Runoff_mm_year;  

    soil_temp     = 10 + 5*sin((t_step-1)/numYears*2*pi);  
    soil_moisture = Soi_mm / max_Soi_mm;                 
    precip         = nanmean(current_precipitation_mm_year(:));
    pet            = nanmean(current_pet_mm_year(:));
    

    [SOC_decomp, DON_min, POP_min, Sorg_min] = ...
        biogeochem_rates(SOC, DON, POP, Sorg, soil_temp, soil_moisture);
    
    [litter_C, litter_N, litter_P, litter_S, uptake_N, uptake_P, uptake_S] = ...
        vegetation_fluxes(SOC, DON, POP, Sorg, soil_temp, soil_moisture, precip, pet);
    

    SOC   = SOC   - SOC_decomp + litter_C;
    DON   = DON   - DON_min    + litter_N;
    POP   = POP   - POP_min    + litter_P;
    Sorg  = Sorg  - Sorg_min   + litter_S;
    
    DIN    = DIN    + DON_min    - uptake_N;
    Pi     = Pi     + POP_min    - uptake_P;
    Sinorg = Sinorg + Sorg_min   - uptake_S;

    SOC(SOC < 0) = 0;
    DON(DON < 0) = 0;
    SIC(SIC < 0) = 0;
    DIN(DIN < 0) = 0;
    POP(POP < 0) = 0;
    Pi(Pi < 0) = 0;
    Sorg(Sorg < 0) = 0;
    Sinorg(Sinorg < 0) = 0;

    runoffMask = Runoff_volume_m3_year > cnpsParams.epsilon;
    
    [DIN, Pi, Sinorg, DIN_leach_ton_year, Pi_leach_ton_year, Sinorg_leach_ton_year] = ...
        cnps_leaching(DIN, Pi, Sinorg, Runoff_volume_m3_year, runoffMask, cnpsParams);

   
    DIN_sim_ton_year(t_step)    = nansum(DIN_leach_ton_year(:));
    Pi_sim_ton_year(t_step)     = nansum(Pi_leach_ton_year(:));
    Sinorg_sim_ton_year(t_step) = nansum(Sinorg_leach_ton_year(:));


    [Cd_ads, Hg_ads, As_ads] = metal_adsorption(SOC_decomp, metalParams);
    [Cd_leach, Hg_leach, As_leach] = metal_leaching(Cd, Hg, As, metalParams);
    [Cd_loss, Hg_loss, As_loss] = metal_erosion_loss( ...
        Cd, Hg, As, ...
        E_water_ha_year, E_wind_ha_year, soil_eroded_underground_ha_year, ...
        metalParams);
    
    Cd = Cd + Cd_precip_input_mg_year + Cd_baseflow_input_mg_year ...
         - Cd_leach - Cd_loss + Cd_ads;
    Hg = Hg + Hg_precip_input_mg_year + Hg_baseflow_input_mg_year ...
         - Hg_leach - Hg_loss + Hg_ads;
    As = As + As_precip_input_mg_year + As_baseflow_input_mg_year ...
         - As_leach - As_loss + As_ads;


    Cd_concentration_3D(t_step,1,:,:) = Cd;
    Hg_concentration_3D(t_step,1,:,:) = Hg;
    As_concentration_3D(t_step,1,:,:) = As;

    concCd = squeeze(Cd_concentration_3D(t_step,:,:,:));  
    concHg = squeeze(Hg_concentration_3D(t_step,:,:,:));
    concAs = squeeze(As_concentration_3D(t_step,:,:,:));

    Cd_concentration_3D(t_step,:,:,:) = metal_diffuse(concCd, diffusion_Cd, delta_z, dt);
    Hg_concentration_3D(t_step,:,:,:) = metal_diffuse(concHg, diffusion_Hg, delta_z, dt);
    As_concentration_3D(t_step,:,:,:) = metal_diffuse(concAs, diffusion_As, delta_z, dt);


    some_threshold = ;  
    pool_matrix = ocean | ((DEM > ZBeachLevel) & (DEM <= ZBeachLevel + some_threshold));


    num_pool_pixels = sum(pool_matrix(:));
    
    CNPS_input_ton_year = DIN_leach_ton_year(pool_matrix) + Pi_leach_ton_year(pool_matrix) + Sinorg_leach_ton_year(pool_matrix); 


    CNPS_input_matrix = zeros(numRows, numCols);
    CNPS_input_matrix(pool_matrix) = CNPS_input_ton_year;  


    CNPS_water_history_ton_year(t_step, :, :) = reshape(CNPS_input_matrix, 1, numRows, numCols); 

 
    Metals_input_mg_year = Cd_leach(pool_matrix) + Hg_leach(pool_matrix) + As_leach(pool_matrix); 


    Metals_input_matrix = zeros(numRows, numCols);
    Metals_input_matrix(pool_matrix) = Metals_input_mg_year; 


    Metals_water_history_mg_year(t_step, :, :) = reshape(Metals_input_matrix, 1, numRows, numCols);  


    water_inflow_ton_year = sum(Runoff_mass_ton_year(:));  
    water_outflow_ton_year = water_inflow_ton_year;    

    total_water_inflow_ton_year = total_water_inflow_ton_year + water_inflow_ton_year;
    total_water_outflow_ton_year = total_water_outflow_ton_year + water_outflow_ton_year;

    oxidation_rate = ; 
    reduction_rate = ;  

    O2_consumption = O2 .* oxidation_rate;  
    O2 = O2 - O2_consumption;


    OR_generation = (O2_consumption .* pixel_area_ha) * (M_O / 1000 / 1000);  
    OR = OR + OR_generation;

    Zmin = min(DEM(:),[],'omitnan');
    Zmax = max(DEM(:),[],'omitnan');
    ZBeachLevel = Zmin + oceanLevelParameter * (Zmax - Zmin);

    ocean = (DEM <= ZBeachLevel) | border_mask;
    VOcean_m3 = nansum((ZBeachLevel - DEM) .* ocean, 'all') * pixel_area; 
    AOcean_m2 = sum(ocean(:)) * pixel_area; 

    [~, linearIdx_new] = min(DEM(:),[],'omitnan');
    [xLow_new, yLow_new] = ind2sub(size(DEM), linearIdx_new);

    if xLow_new ~= xLow || yLow_new ~= yLow
        drain(:) = false;
        drain(xLow_new, yLow_new) = true;
        xLow = xLow_new;
        yLow = yLow_new;
    end

    AET_sim(t_step) = nanmean(current_pet_mm_year(:));  

    if t_step == 1
        Runoff_sim = zeros(numRows, numCols, numYears, 2);  
        Baseflow = zeros(numRows, numCols, numYears);      
    end

    Runoff_sim(:,:,t_step,1) = Runoff_mm_year;        
    Runoff_sim(:,:,t_step,2) = Runoff_mass_ton_year;   
    Baseflow(:,:,t_step) = current_baseflow_kg_m2_year; 

    E_water_history_ha_year(t_step, :, :) = E_water_ha_year;        
    E_wind_history_ha_year(t_step, :, :) = E_wind_ha_year;          
    E_total_history_ha_year(t_step, :, :) = total_erosion_ha_year;      

    D_history_ha_year(t_step, :, :) = D_ha_year;                    
    SOC_history_ha_year(t_step, :, :) = SOC;          
    SIC_history_ha_year(t_step, :, :) = SIC;         
    DON_history_ha_year(t_step, :, :) = DON;          
    DIN_history_ha_year(t_step, :, :) = DIN;          
    POP_history_ha_year(t_step, :, :) = POP; 

    Pi_history_ha_year(t_step, :, :) = Pi;           
    Sorg_history_ha_year(t_step, :, :) = Sorg;        
    Sinorg_history_ha_year(t_step, :, :) = Sinorg;    

    O2_history_mg_L_year(t_step, :, :) = O2;            
    OR_history_mol_O_ha_year(t_step, :, :) = OR;       
    Cd_history_mg_ha_year(t_step, :, :) = Cd;         
    Hg_history_mg_ha_year(t_step, :, :) = Hg;          
    As_history_mg_ha_year(t_step, :, :) = As;          

    soil_generated_metals_history_mg_ha_year(t_step, :, :, 1) = soil_generated_Cd; 
    soil_generated_metals_history_mg_ha_year(t_step, :, :, 2) = soil_generated_Hg;  
    soil_generated_metals_history_mg_ha_year(t_step, :, :, 3) = soil_generated_As; 
    
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 1) = soil_eroded_Cd;     
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 2) = soil_eroded_Hg;     
    soil_eroded_metals_history_mg_ha_year(t_step, :, :, 3) = soil_eroded_As;     


    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 1) = E_water_ha_year .* metal_eroded_ratio_Cd * 1000;  
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 2) = E_water_ha_year .* metal_eroded_ratio_Hg * 1000; 
    soil_eroded_water_metals_history_mg_ha_year(t_step, :, :, 3) = E_water_ha_year .* metal_eroded_ratio_As * 1000; 

    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 1) = E_wind_ha_year .* metal_eroded_ratio_Cd * 1000;   
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 2) = E_wind_ha_year .* metal_eroded_ratio_Hg * 1000;  
    soil_eroded_wind_metals_history_mg_ha_year(t_step, :, :, 3) = E_wind_ha_year .* metal_eroded_ratio_As * 1000;  
end
