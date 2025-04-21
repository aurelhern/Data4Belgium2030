include("./input_MC/sampling_functions_v6.jl")
include("./input_MC/practical_functions.jl")

function MC(simu_name, excel_file_input_data, n_min, n_max, n_threads, n_batch, n_years_in_batch, print_yearlydata, print_verifyelec, profile, import_avail, timestep, outagestech, outageswind, constance, incr_limit, convergnce_type, seed, min_lf_el)

    start_simu = now()

    seed1 = seed[1]
    seed2 = seed[2]
    seed3 = seed[3]
    seed4 = seed[4]
    seed5 = seed[5]
    seed6 = seed[6]

## 0. Save input data
    file_output = "./Results/"*simu_name*"/output/"
    xlsx_file_data = XLSX.open_xlsx_template(excel_file_input_data)
    XLSX.writexlsx(file_output*"input_data.xlsx", xlsx_file_data, overwrite=true)

## 1. Set up Static Parameters ##

    ## Read Grid Topology Parameters ##

    Data = XLSX.readxlsx(excel_file_input_data)

    Buses = DataFrame(XLSX.gettable(Data["Buses"]))
    Gen = DataFrame(XLSX.gettable(Data["Gen"]))
    Lines = DataFrame(XLSX.gettable(Data["Lines"]))
    Res = DataFrame(XLSX.gettable(Data["Res"]))
    H2PipeLines = DataFrame(XLSX.gettable(Data["H2PipeLines"]))
    Electrolyzer = DataFrame(XLSX.gettable(Data["Electrolyzer"]))
    H2toP = DataFrame(XLSX.gettable(Data["H2toP"]))
    H2Sto = DataFrame(XLSX.gettable(Data["H2Sto"]))
    Imp = DataFrame(XLSX.gettable(Data["Imp"]))
    Phs = DataFrame(XLSX.gettable(Data["Phs"]))
    Batt = DataFrame(XLSX.gettable(Data["Batt"]))
    General = DataFrame(XLSX.gettable(Data["General"]))
    H2Imports = DataFrame(XLSX.gettable(Data["H2Imports"])) 

    res_nidx = indexin(Res[!, :Node], Buses[!, :Name])
    phs_nidx = indexin(Phs[!, :Node], Buses[!, :Name])
    batt_nidx = indexin(Batt[!, :Node], Buses[!, :Name])
    imp_nidx = indexin(Imp[!, :Node], Buses[!, :Name])
    g_nidx = indexin(Gen[!, :Node], Buses[!, :Name])
    l_fromnidx = indexin(Lines[!, :FromNode], Buses[!, :Name])
    l_tonidx = indexin(Lines[!, :ToNode], Buses[!, :Name])

    h2pipe_fromnidx = indexin(H2PipeLines[!, :FromNode], Buses[!, :Name])
    h2pipe_tonidx = indexin(H2PipeLines[!, :ToNode], Buses[!, :Name])
    el_nidx = indexin(Electrolyzer[!, :Node], Buses[!, :Name])
    h2p_nidx = indexin(H2toP[!, :Node], Buses[!, :Name])
    h2sto_nidx = indexin(H2Sto[!, :Node], Buses[!, :Name])
    h2imp_nidx = indexin(H2Imports[!, :Node], Buses[!, :Name])

    resn = Dict{Int, Vector{Int}}() #for each bus n : list of renewable assets connected to it 
    phsn = Dict{Int, Vector{Int}}() #for each bus n : list of phs connected to it 
    battn = Dict{Int, Vector{Int}}() #for each bus n : list of batteries connected to it 
    impn = Dict{Int, Vector{Int}}() #for each bus n : list of international interconnections connected to it 
    gn = Dict{Int, Vector{Int}}() #for each bus n : list of generators connected to it 
    l_fromn = Dict{Int, Vector{Int}}() #for each bus n : list of (from) lines 
    l_ton = Dict{Int, Vector{Int}}() #for each bus n : list of (to) lines 
    h2pipe_fromn = Dict{Int, Vector{Int}}() #for each bus n : list of (from) h2 pipelines 
    h2pipe_ton = Dict{Int, Vector{Int}}() #for each bus n : list of (to) h2 pipelines 
    eln = Dict{Int, Vector{Int}}() #for each bus n : list of electrolyzers connected to it 
    h2pn = Dict{Int, Vector{Int}}() #for each bus n : list of H2P unit connected to it 
    h2ston = Dict{Int, Vector{Int}}() #for each bus n : list of H2 storage unit connected to it 
    h2impn = Dict{Int, Vector{Int}}() #for each bus n : list of H2 import connections connected to it 


    for n = 1:size(Buses, 1)
        resn[n] = Int[]
        phsn[n] = Int[]
        battn[n] = Int[]
        impn[n] = Int[]
        gn[n] = Int[]
        l_fromn[n] = Int[]
        l_ton[n] = Int[]
        h2pipe_fromn[n] = Int[]
        h2pipe_ton[n] = Int[]
        eln[n] = Int[]
        h2pn[n] = Int[]
        h2ston[n] = Int[]
        h2impn[n] = Int[]
    end
    for r = 1:size(Res, 1)
        push!(resn[res_nidx[r]], r)
    end
    for p = 1:size(Phs, 1)
        push!(phsn[phs_nidx[p]], p)
    end
    for b = 1:size(Batt, 1)
        push!(battn[batt_nidx[b]], b)
    end
    for g = 1:size(Gen, 1)
        push!(gn[g_nidx[g]], g)
    end
    for i = 1:size(Imp, 1)
        push!(impn[imp_nidx[i]], i)
    end
    for l = 1:size(Lines, 1)
        push!(l_fromn[l_fromnidx[l]], l)
        push!(l_ton[l_tonidx[l]], l)
    end
    for l = 1:size(H2PipeLines, 1)
        push!(h2pipe_fromn[h2pipe_fromnidx[l]], l)
        push!(h2pipe_ton[h2pipe_tonidx[l]], l)
    end
    for e = 1:size(Electrolyzer, 1)
        push!(eln[el_nidx[e]], e)
    end
    for h = 1:size(H2toP, 1)
        push!(h2pn[h2p_nidx[h]], h)
    end
    for j = 1:size(H2Sto, 1)
        push!(h2ston[h2sto_nidx[j]], j)
    end
    for i = 1:size(H2Imports, 1)
        push!(h2impn[h2imp_nidx[i]], i)
    end

    ## Read Techno-Economic Data of technologies ##

    voll= reshape(General[!,"VOLL_€/MWh"],1)[1] # Value of Loss-Load [€/MWh] 
    voll_h2 = voll/10 # Value of Loss-Load [€/MWh] 
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] ;# Days per months [days/months]
    sb = reshape(General[!,"SB_MVA"],1)[1] # Base power of the per-unit system [MVA]
    nref = reshape(General[!,"RefBus"],1)[1] # Refernce node (voltage angle Theta = 0)
    lf_chp_min = reshape(General[!,"LfChp"],1)[1] # Constant load factor of individual CHP units [-]

    # Conventional Generators

    pconv_max = Gen[!, :Pmax_MW]  # Max capacity of each generator [MW]

    Gen[occursin.("CHP", Gen.Type), :Pmin_MW] = lf_chp_min*Gen[occursin.("CHP", Gen.Type), :Pmax_MW] # Min capacity of individual CHP units [MW]
    Gen[occursin.("Biomass", Gen.Fuel), :Pmin_MW] = Gen[occursin.("Biomass", Gen.Fuel), :Pmax_MW] # Min capacity of individual biomass units [MW]
    Gen[occursin.("Waste", Gen.Fuel), :Pmin_MW] = Gen[occursin.("Waste", Gen.Fuel), :Pmax_MW] # Min capacity of individual waste incinerators units [MW]
    Gen[occursin.("NU", Gen.Type), :Pmin_MW] = Gen[occursin.("NU", Gen.Type), :Pmax_MW] # Min capacity of nuclear powr plants [MW]
    idx_nu = Gen[occursin.("NU", Gen.Type), :Generator] # Index of nuclear units

    pconv_min = Gen[!, :Pmin_MW];  # Min capacity of each generator [MW]

    mttf_conv = Gen[!, :MTTF] # MTTF of each generator [hours/year]
    mttr_conv = Gen[!, :MTTR] # MTTFRof each generator [hours/year]
    
    c_emissions = reshape(General[!,"CarbonCost_€/tonCO2"],1)[1] # Carbon Cost [€/tonCO2]
    c_gen_fuel = Gen[!,"FuelCost_€/MWh_fuel"] # Fuel cost [€/MWh_fuel]
    emission_carbon_fuel = Gen[!,"CO2Emissions_tCO2/MWh_fuel"] # CO2 Emissions [tCO2/MWh_fuel]
    bool_taxcarbon = [parse(Bool, x) for x in Gen[!,"TaxCO2Emissions"]]  # Unit concerned by CO2 tax? [boolean]
    c_gen_vom = Gen[!,"VOM_€/MWh"] # VOM [€/MWh]
    eff_gen = Gen[!,"Efficiency"] # Generator Efficiency [-]
    c_gen_op = c_gen_fuel./eff_gen + bool_taxcarbon.*emission_carbon_fuel./eff_gen.*c_emissions + c_gen_vom # Marginal cost of each generator [€/MWh]

    open("./Results/"*simu_name*"/costgenop.txt", "w") do io
        writedlm(io, round.(c_gen_op, digits=4))   # print operation cost of each generator
    end

    # Interconnections
    pimp_max = Imp[!, "Pmax_MW"] # Import capacity of each interconnection [MW]
    eimp_year_max = 1e6*Imp[!, "Emax_year_TWh"] # Annual import capacity of each interconnection [MWh/year] 
    pimp_max_tot = reshape(General[!,"MaxPowerImports_MW"],1)[1] #Simultaneous import limit of electricity [MW]

    import_countries = Imp[!, "Country"] # Country of importation
    if size(Imp)[1] > 0 # If interconnections
        if import_avail == "random"
            c_imp = Imp[!, "Price_€/MWh"] # Price of imports at each interconnection [€/MWh]
        elseif import_avail == "max"          
            c_imp_base = import_price(import_countries, timestep)
            mean_noise_c_imp = Imp[!, "Mean_noise_€/MWh"]
            std_noise_c_imp = Imp[!, "Std_noise_€/MWh"]
            seed_noise_c_imp = Int.(Imp[!, "Seed_noise"])
        end   
    else
        c_imp = 0
    end

    mttf_imp = Imp[!, "MTTF"] # MTTF of each interconnection [hours/year]
    mttr_imp = Imp[!, "MTTR"]; # MTTR of each interconnection[hours/year]

    # Pumped Hydro Storage (PHS)
    eff_pump = Phs[!, "eff_pump"] # Efficiency of PHS pumps [-] 
    eff_turb = Phs[!, "eff_turb"] # Efficiency of PHS turbines [-] 
    soc_phs_max = Phs[!, "Emax_MWh"] # PHS storage capacity [MWh_el] 
    p_phs_max = Phs[!, "Pmax_MW"] #Power of PHS  pumps/turbines [MW] 

    # Batteries
    eff_batt_in = Batt[!, "eff_batt_in"] # Input efficiency of battery 
    eff_batt_out = Batt[!, "eff_batt_out"] # Output efficiency of battery 
    soc_batt_max = Batt[!, "Emax_MWh"]  # Battery storage capacity [MWh_el] 
    p_batt_max = Batt[!, "Pmax_MW"] # Input/Output power of batteries 

    # Renewable Energy Sources
    if size(Res)[1] > 0 # If offshore wind farms connected to the HV grid
        pwindoffshore_max = sum(Res[!, :Pmax_MW]) #Total power of offshore wind farms [MW] 
    else # If no renewable power plant connected to the HV grid
        pwindoffshore_max = 0 
    end

    wf_names = Res[!, "Name"] # Name of offshore wind farms
    nb_wt_offshore = Res[!, "Nb_wt"] # Number of wind turbines in each offshore wind farm
    mttf_res = Res[!, "MTTF"] # MTTF of wind turbines in the concerned wind farm
    mttr_res = Res[!, "MTTR"] # MTTR of wind turbines in the concerned wind farm
    increasefactoroffshore = pwindoffshore_max/reshape(General[!,"PowerWfThnCode_MW"],1)[1] # Increase factor from wind farm power generator of THN 
    println("Increase Factor: ", increasefactoroffshore)

    #Extract voltages of Buses
    ub_l_from = zeros(size(Lines,1))
    ub_l_to = zeros(size(Lines,1))
    ub_l = zeros(size(Lines,1))
    for l = 1:size(Lines,1)
        ub_l_from[l] = Buses[l_fromnidx[l],"UB_kV"]
        ub_l_to[l] = Buses[l_tonidx[l],"UB_kV"]
        ub_l[l] = max(ub_l_from[l], ub_l_to[l])
    end

    # Lines 
    f_max = Lines[!, "Pmax_MW"] # Power capacity of the line [MW]
    angmin_rad = -30/360*2*pi # Minimum voltage angle of bus [rad]
    angmax_rad = 30/360*2*pi # Maximum voltage angle of bus [rad]
    mttf_f = Lines[!, "MTTF"] # MTTF of each line [hours/year]
    mttr_f = Lines[!, "MTTR"] # MTTR of each line [hours/year]
    M = 100000000000000 #infinite value

    # Lines susceptance
    x_ohm = Lines[!,"X_ohm"]
    zb = ub_l.^2/sb # Base impedence of HV grid [ohm]    
    x_pu = x_ohm./zb # Base impedence of HV grid [pu]    
    b_s = x_ohm.^(-1) # Susceptance of lines [S]
    b_pu = x_pu.^(-1) # Susceptance of lines [pu]

    # Space repartition 
    d_l = 0.01*Buses[!,"Load_share_%"] # Load share at each bus [-]
    d_pv = 0.01*Buses[!,"Pv_share_%"] # PV share at each bus [-]
    d_onshore = 0.01*Buses[!,"Onshore_share_%"] # Onshore wind share at each bus [-]
    d_hydro = 0.01*Buses[!,"Hydro_share_%"] # Hydro share at each bus [-]
    d_hp = 0.01*Buses[!,"Hp_share_%"] # Heat-pumps share at each bus [-]
    d_ev = 0.01*Buses[!,"Ev_share_%"] # Electric vehicles share at each bus [-]
    d_sh = 0.01*Buses[!,"Load_share_%"] # Sheddable load share at each bus [-]
    d_s = 0.01*Buses[!,"Load_share_%"] # Shiftable load share at each bus [-]
    d_h2 = 0.01*Buses[!,"H2_share_%"] # H2 load share at each bus [-]
    d_distgen = 0.01*Buses[!,"DistGen_share_%"] # DistGen share at each bus [-]

    # Electric Vehicles 
    EVs = DataFrame(XLSX.gettable(Data["EVs"])) #Dataset of EVs

    n_ev = reshape(EVs[!,"NumberEvs"],1)[1] # Number of EVs in 2030 [cars]

    ev_v0 = readdlm("./input_MC/V0.txt")/100  #Daily distribution [%/h/car] - vector of 24 values
    ev_v1h = readdlm("./input_MC/V1H.txt")/100 #Daily distribution [%/h/car] - vector of 24 values
    ev_v2h = readdlm("./input_MC/V2H.txt")/100 #Daily distribution [%/h/car] - vector of 24 values

    ev_v0 = ev_v0/sum(ev_v0)
    ev_v1h = ev_v1h/sum(ev_v1h)
    ev_v2h = ev_v2h/sum(ev_v2h)

    w_v0 = 0.01*reshape(EVs[!,"V0_share_%"],1)[1] # Share of V0 vehicles [-]
    w_v1h = 0.01*reshape(EVs[!,"V1H_share_%"],1)[1] # Share of V1H vehicles [-]
    w_v2h = 0.01*reshape(EVs[!,"V2H_share_%"],1)[1] # Share of V2H vehicles [-]

    w_v1m = 0.01*reshape(EVs[!,"V1M_share_%"],1)[1] # Share of V2H vehicles [-]
    avail_v1m = repeat(readdlm("./input_MC/avail_V1M.txt")/100,365) #Daily availability V1M [%/h] - vector of 24 values
    p_ev_charger = reshape(EVs[!,"ChargerPower_kW"],1)[1]/1000 # Power of EV charger [MW]

    yearly_km = reshape(EVs[!,"YearlyDistance_km/year/car"],1)[1] # [km/year/car]
    eff_car = reshape(EVs[!,"CarEfficiency_kWh/km"],1)[1] # [kwh/km]
    yearly_need = yearly_km*eff_car # [kwh/year/car]
    daily_need = yearly_need/365 # [kwh/day/car]

    ev_charging_daily = n_ev*(w_v0*ev_v0 + w_v1h*ev_v1h + w_v2h*ev_v2h)*daily_need/1000 # Resulting daily weighted charging profile [MWh/h] for 24 hours - vector of 24 values
    load_ev_fixed = repeat(ev_charging_daily,365) # Fixed EV Load [MWh/h] for 8760 hours

    load_ev_flex_daily = n_ev*w_v1m*daily_need/1000 # Daily Flex EV Load [MWh/day] 

    # Heat Pumps
    HPs = DataFrame(XLSX.gettable(Data["HPs"])) #Dataset of EVs

    hp0 = readdlm("./input_MC/HP0.txt")  # Daily distribution [%/h/car]- vector of 24 values
    hp1h = readdlm("./input_MC/HP1H.txt") # Daily distribution [%/h/car]- vector of 24 values

    hp0 = hp0/sum(hp0)
    hp1h = hp1h/sum(hp1h)

    w_hp0 = 0.01*reshape(HPs[!,"HP0_share_%"],1)[1] # Share of HP0 heat-pumps [-]
    w_hp1h = 0.01*reshape(HPs[!,"HP1H_share_%"],1)[1] # Share of HP1H heat-pumps [-]
    w_hp1m = 0.01*reshape(HPs[!,"HP1M_share_%"],1)[1] # Share of HP1M heat-pumps [-]

    hp_monthly_dist = readdlm("./input_MC/hp_monthly_dist.txt") # Daily distribution [%/h/car] - vector of 12 values
    hp_yearly_demand = reshape(HPs[!,"YearlyHpDemand_TWh/year"],1)[1]*1e6 # Yearly electrical demand from heat-pumps [MWh/year]

    load_hp_fixed = ones(8760)
    for m=1:12
        if m==1
            start = 1
        end
        if m>1
            start =  sum(days_per_month[i] for i=1:m-1)*24+1
        end
        stop = sum(days_per_month[i] for i=1:m)*24
        load_hp_fixed[start:stop] = hp_monthly_dist[m]*hp_yearly_demand/days_per_month[m]*repeat(w_hp0*hp0 + w_hp1h*hp1h, days_per_month[m]) # Hourly electrical demand of heat pumps [MWh_el/h] - vector of 8760 values
    end

    load_hp_flex = w_hp1m*[hp_monthly_dist[m]*hp_yearly_demand/days_per_month[m] for m=1:12] # Daily demand of HPs for HP1M [MWh/day] - vector of 12 values
    a = [load_hp_flex[m]*ones(days_per_month[m]) for m=1:12]
    load_hp_flex = vcat(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12]) # Daily demand of HPs for HP1M [MWh/day] - vector of 365 values

    n_hp = reshape(HPs[!,"NumberHps"],1)[1] # Number of HPs in 2030 [hps]
    cop_hp = reshape(HPs[!,"COP"],1)[1] # COP of heat-pumps [W_th/W_el]
    p_hp_individual = reshape(HPs[!,"Pmax_kWth"],1)[1]/1000/cop_hp # Elec power of installed HP [MW_el]

    open("./Results/"*simu_name*"/load_hp_flex.txt", "w") do io
        writedlm(io, round.(load_hp_flex, digits=4))   # print load_hp_flex
    end

    open("./Results/"*simu_name*"/load_hp_fixed.txt", "w") do io
        writedlm(io, round.(load_hp_fixed, digits=4))   # print load_hp_fixed
    end

    #Distributed Renewable Power
    DistRes = DataFrame(XLSX.gettable(Data["DistRes"])) #Dataset of heat consumption profiles 

    psolar = reshape(DistRes[!,"Psolar_MW"],1)[1] # PV capacity [MW] 
    ponshore = reshape(DistRes[!,"Ponshore_MW"],1)[1] # Wind onshore capacity [MW] 
    phydro = reshape(DistRes[!,"Phydro_MW"],1)[1] # Hydro capacity [MW] 

    solar_lf = readdlm("./input_MC/data_distres/Time_series_PV.txt")[1:timestep] # PV load factor [-] - vector of 8760 values
    onshore_lf = readdlm("./input_MC/data_distres/Time_series_Wind_onshore.txt")[1:timestep] # Onshore wind power factor [-] - vector of 8760 values
    hydro_lf = readdlm("./input_MC/data_distres/Time_series_Hydro_river.txt")[1:timestep] # Onshore wind power factor [-] - vector of 8760 values

    #Elec load
    yearly_demand_elec = reshape(General[!,"YearlyDemand_TWh/year"], 1)[1]*1e6 #Yearly electrical demand [MWh/year]
    
    load_elia_pattern = readdlm("./input_MC/output_data_process/eliagridload_mw_2018.txt") 
    load_elia_pattern = load_elia_pattern/sum(load_elia_pattern)
    
    load = load_elia_pattern*yearly_demand_elec

    #Distributed generators
    DistGen = DataFrame(XLSX.gettable(Data["DistGen"]))

    p_distgen_max_chp = DistGen[occursin.("Gas CHP", DistGen.Type), :Pmax_MW] # Aggregated power of distributed gas CHPs [MW_el]
    p_distgen_max_biomass = DistGen[occursin.("Biomass", DistGen.Type), :Pmax_MW] # Aggregated power of distributed biomass plants [MW_el]
    p_distgen_max_waste = DistGen[occursin.("Waste", DistGen.Type), :Pmax_MW] # Aggregated power of distributed waste incinerators [MW_el]

    lf_distgen_biomass = ones(8760).*DistGen[occursin.("Biomass", DistGen.Type), :LF] # Load factor of distributed gas CHPs [-]
    lf_distgen_waste = ones(8760).*DistGen[occursin.("Waste", DistGen.Type), :LF] # Load factor of distributd biomass plants [-]
    lf_distgen_chp = load/sum(load)*8760*DistGen[occursin.("Gas CHP", DistGen.Type), :LF] # Load factor of distributed waste plants [-]

    pdistgen = p_distgen_max_chp.*lf_distgen_chp + p_distgen_max_waste.*lf_distgen_waste + p_distgen_max_biomass.*lf_distgen_biomass # Aggregated operation power of distributed generation [MW_el]

    # Hydrogen Demand
    H2demand = DataFrame(XLSX.gettable(Data["H2demand"])) #Dataset of hydrogen demand

    lhv_H2 = 33.33 # H2 LHV [MWh_H2/t_H2]
    qh2_demand_year_kt = reshape(H2demand[!, "YearlyH2Demand_kt/year"],1)[1] # Yearly H2 demand [kt_H2/year]
    qh2_demand_year_mwh = qh2_demand_year_kt*1000*lhv_H2 # Yearly H2 demand [MWh_H2/year]

    w_h2flex = 0.01*reshape(H2demand[!,"H2Flex_share_%"],1)[1] # Share of flexible H2 demand [-]
    w_h2fixed = 0.01*reshape(H2demand[!,"H2Fixed_share_%"],1)[1] # Share of fixed H2 demand [-]

    load_h2_fixed = w_h2fixed*qh2_demand_year_mwh/8760*ones(timestep) # Yearly H2 demand distributed uniformly throughout the year [MWh_H2/h] - vector of 8760 values
    load_h2_flex = w_h2flex*qh2_demand_year_mwh # Yearly flexible H2 demand [MWh_H2/year]

    # H2 Technologies

    ph2_electrolyzer_max = Electrolyzer[!, :Pmax_MW] # Installed power of electrolyzers [MW_el] 
    mttf_el = Electrolyzer[!, :MTTF] # MTTF of each electrolyzer [hours/year]
    mttr_el = Electrolyzer[!, :MTTR] # MTTR of each electrolyzer [hours/year]
    eff_electrolyzer =  0.7 #Electrolyzer[!, "Efficiency"] # Efficiency of electrolyzers  [-]
    compr_H2 = 0.09 #Electrolyzer[!, "CompressionEfficiency"] # Electrical consumption for H2 compression to 350 bars [MWh_el/MWh_H2] 
    cel = Electrolyzer[!, "CostElectrolyzer_€/MWh"] # Operational cost of electrolyzers (not including the fuel cost and compression needs) [€/Mwh_el]

    
    p_h2top_max = H2toP[!, :Pmax_MW] # Installed power of hydrogen-to-power units [MW_el] 
    mttf_h2p = H2toP[!, :MTTF] # MTTF of each hydrogn-to-power unit [hours/year]
    mttr_h2p = H2toP[!, :MTTR] # MTTR of each hydrogn-to-power unit [hours/year]
    eff_h2top = 0.43 #H2toP[!, "Efficiency"] # Efficiency of hydrogen-to-power units [-]

    soc_h2_max = H2Sto[!,"Emax_MWh_H2"] # Maximum capacity of H2 storage unit [MWh_h2]
    q_h2sto_in_max = H2Sto[!,"Pmax_MW"] # Maximum input flows of H2 storage unit [MWh_h2/h]
    q_h2sto_out_max = H2Sto[!,"Pmax_MW"] # Maximum output flows of H2 storage unit [MWh_h2/h]
    eff_h2sto_in = H2Sto[!,"eff_h2sto_in"] # Input Efficiency of H2 Storage [-]
    eff_h2sto_out = H2Sto[!,"eff_h2sto_out"] # Output Efficiency of H2 Storage [-]

    f_h2_max = H2PipeLines[!,"Qmax_MWh_h2/h"] # Maximum h2 flow through pipelines [MWh_h2/h]

    qh2_import_max = H2Imports[!,"Qmax_MWh_h2/h"]  # Maximum h2 import flow [MWh_h2/h]
    ch2imp = H2Imports[!,"PriceH2Imports_€/MWh"]  # Price of H2 imports [€/Mwh_h2]

    # Load Shifting
    Shift = DataFrame(XLSX.gettable(Data["Shift"])) # Dataset of load shifting 

    premove_max = reshape(Shift[!,"Pmax_MW"],1)[1] # Max shiftable power [MW_el]
    sremove_day_max = reshape(Shift[!,"Emax_daily_MWh/day"],1)[1] # Maximum daily shiftable energy [MWh_el/day]
    cshift = reshape(Shift[!,"Price_€/MWh"],1)[1] # price to shift energy [€/Mwh]

    # Voluntary Load Sheding
    Shed = DataFrame(XLSX.gettable(Data["Shed"])) #Dataset of load sheding 

    pshed_max = reshape(Shed[!,"Pmax_MW"],1)[1] # Max voluntary shedable power [MW]
    eshed_day_max = reshape(Shed[!,"Emax_daily_MWh/day"],1)[1] # Max voluntary shed energy per day [MWh/day]
    cshed = reshape(Shed[!,"Price_€/MWh"],1)[1] # Price to shed energy [€/MWh]

## 2. Initial Dynamic Parameters ##

    # Availabilities due to potential outages
    if outagestech==true # If outages are taken into account 

        availconv = avail_comp(length(pconv_max), timestep, mttf_conv, mttr_conv, pconv_max, false, zeros(timestep), seed1)[1] # Availability of Generators (TTF and TTR sampling)
        availline = avail_comp(lastindex(f_max), timestep, mttf_f, mttr_f, f_max, false, zeros(timestep), seed2)[1] # Availability of Lines (TTF and TTR sampling)

        if size(Electrolyzer)[1] > 0
            availel = avail_comp(lastindex(ph2_electrolyzer_max), timestep, mttf_el, mttr_el, ph2_electrolyzer_max, false, zeros(timestep), seed3)[1] # Availability of electrolyzers (TTF and TTR sampling)
        end

        if size(H2toP)[1] > 0
            availh2p = avail_comp(lastindex(p_h2top_max), timestep, mttf_h2p, mttr_h2p, p_h2top_max, false, zeros(timestep), seed4)[1] # Availability of hydrogen-to-power units (TTF and TTR sampling)
        end
    else
        availconv = ones(length(pconv_max), timestep) # Availability of Generators
        availline = ones(lastindex(f_max), timestep) # Availability of Lines
        availel = ones(lastindex(ph2_electrolyzer_max), timestep) # Availability of electrolyzers
        availh2p = ones(lastindex(p_h2top_max), timestep) # Availability of hydrogen-to-power units
    end

    # Availabilities of interconnections and neighbouring generation potential
    if size(Imp)[1] > 0 # If interconnections
        if import_avail == "random"
            availimp = 2*rand(size(Imp,1), timestep)
            replace(x -> x>1 ? 0 : x, availimp)
        elseif import_avail == "max"
            availimp = avail_comp(lastindex(pimp_max), timestep, mttf_imp, mttr_imp, pimp_max, false, zeros(timestep), seed5)[1]
            c_imp = Array{Float64}(undef, size(import_countries)[1], timestep)
            for c in 1:size(import_countries)[1]
                Random.seed!(seed_noise_c_imp[c]) # Setting the seed
                dist_c_imp = Normal(mean_noise_c_imp[c], std_noise_c_imp[c])
                noise = rand(dist_c_imp, (timestep))
                c_imp[c, :] = c_imp_base[c, :] + noise
            end
        end   
    else
        c_imp = 0
        availimp = zeros(size(Imp,1), timestep)
    end

    # Offshore wind farms production
    if size(Res)[1] > 0 # If no interconnection
        res_f = offshore_time_series4(nb_wt_offshore, wf_names, timestep, mttf_res, mttr_res, outageswind, constance, seed6) # Offshore wind farm power [MW_el] - matrix: 8760 values for each wind farm
        res = increasefactoroffshore*res_f[1]
        r_wind = res_f[2]
    else
        res = zeros(2,timestep)
    end

## Initialization of Monte-Carlo results

    global epsilon_loee = Vector{Float64}()
    #global epsilon_loee = append!(epsilon_loee,1)
    global epsilon_loee_old = Vector{Float64}()
    global incr_epsilon_loee = Vector{Float64}()
    global incr_epsilon_loee_rm = Vector{Float64}()
    global n_mc=0

    global lole =  Vector{Float64}()

    global lol_set = Vector{Int64}()
    global lol_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)
    global loe_set = Vector{Float64}()
    global loe_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)
    global h2loe_set = Vector{Float64}()
    global h2loe_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)

    global e_curt_set = Vector{Float64}()
    global e_curt_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)
    global p_curt_set = Vector{Float64}()

    global phscycle_set = Array{Float64}(undef, size(Phs)[1], 0)
    global battcycle_set = Array{Float64}(undef, size(Batt)[1], 0)
    global h2cycle_set = Array{Float64}(undef, size(H2Sto)[1], 0)

    global shift_set = Vector{Float64}()
    global shift_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)
    global shed_set = Vector{Float64}()
    global shed_nodal_set = Array{Float64}(undef, size(Buses)[1], 0)

    global conv_set = Vector{Float64}()
    global genlf_set = Array{Float64}(undef, size(Gen)[1], 0)
    global offshore_set = Vector{Float64}()

    global imp_set = Vector{Float64}()
    global imp_nodal_set = Array{Float64}(undef, size(Imp)[1], 0)
    global h2imp_set = Vector{Float64}() 
    global h2imp_nodal_set = Array{Float64}(undef, size(H2Imports)[1], 0)

    global lly_set = Array{Float64}(undef, size(Lines)[1], 0)
    global plly_set = Array{Float64}(undef, size(H2PipeLines)[1], 0)

    global elinst_set = Vector{Float64}()
    global elinst_nodal_set = Array{Float64}(undef, size(Electrolyzer)[1], 0)
    global ellf_set = Vector{Float64}()
    global ellf_nodal_set = Array{Float64}(undef, size(Electrolyzer)[1], 0)

    global h2topinst_set = Vector{Float64}() 
    global h2topinst_nodal_set = Array{Float64}(undef, size(H2toP)[1], 0)
    global h2toplf_set = Vector{Float64}() 
    global h2toplf_nodal_set = Array{Float64}(undef, size(H2toP)[1], 0)

    global soc_h2_max_res_set = Vector{Float64}() 

    global cost_set = Vector{Float64}()
    global runtimeopti_set = Vector{Float64}()

    #Prints parameters
    if print_yearlydata == false
        i_yearlydata = 3
    else
        i_yearlydata = n_max+1
    end

    if print_verifyelec == false
        i_verifyelec = 1
    else
        i_verifyelec = n_max+1
    end

## Thread loop

    u = ReentrantLock()
    @threads for thread in 1:n_threads

## 3. Optimization Model ##

        start_modelopt = now()
        println("--- Thread ", Threads.threadid(),"- Time : "  , hour(start_modelopt), ":", minute(start_modelopt), ":", second(start_modelopt),  " - Start building model") # ~13.8 s

        model = Model(Gurobi.Optimizer)
        # The barrier solver terminates when the relative difference between the primal and dual objective 
        #values is less than the specified tolerance (with a GRB_OPTIMAL status). Tightening this tolerance 
        #often produces a more accurate solution, which can sometimes reduce the time spent in crossover. 
        #Loosening it causes the barrier algorithm to terminate with a less accurate solution, which can 
        #be useful when barrier is making very slow progress in later iterations.
        #MOI.set(model, MOI.RawParameter("FeasibilityTol"), 1e-3) 
        #MOI.set(model, MOI.RawParameter("OptimalityTol"), 1e-3)
        #MOI.set(model, MOI.RawParameter("Crossover"), 0) #Removes the crossover which makes sure the interior solution founf by the barrier method is llocated on a vertex.
        #MOI.set(model, MOI.RawParameter("BarConvTol"), 1e-10) #Default value: 1e-8. The smaller it is, the shorter the crossover will take.

        # if Threads.threadid() != 1
        #     set_optimizer_attributes(model, "OutputFlag" => 0)
        # end
        set_optimizer_attributes(model, "OutputFlag" => 0)
        set_optimizer_attribute(model,"Threads", 5);

        ## Variables

        @variable(model, Pconv[g=1:size(Gen, 1), t=1:timestep] >=0) # Operating power of each generator [MW]
        @variable(model, F[k=1:size(Lines, 1), t=1:timestep]) # Power flowing through each lines [MW]
        @variable(model, Theta[n=1:size(Buses, 1), t=1:timestep]) # Phase angle of each node [rad]
        @variable(model, Ens[n=1:size(Buses,1), t=1:timestep] >=0) # Energy not served at each node [MWh]
        @variable(model, Ens_h2[n=1:size(Buses,1), t=1:timestep] >=0) # H2 Energy not served at each node [MWh]
        @variable(model, Pimp[i=1:size(Imp,1), 1:timestep]>=0) # Enegry imported at each interconnection [MW]
        @variable(model, Curtail[n=1:size(Buses,1), t=1:timestep] >=0) # Curtailed energy at each node [MW]

        @variable(model, Soc_phs[p=1:size(Phs,1), 1:timestep+1] >=0) # State-of-charge of PHS [MWh_el]
        @variable(model, Pphs_sto_in[p=1:size(Phs,1),1:timestep+1] >=0) # Input power of PHS storage [MW_el]
        @variable(model, Pphs_sto_out[p=1:size(Phs,1),1:timestep+1] >=0) # Output power of PHS storage [MW_el]

        @variable(model, Soc_batt[p=1:size(Batt,1), 1:timestep+1] >=0) # State-of-charge of batteries [MWh_el]
        @variable(model, Pbatt_sto_in[p=1:size(Batt,1),1:timestep+1] >=0) # Input power of batteries [MW_el]
        @variable(model, Pbatt_sto_out[p=1:size(Batt,1),1:timestep+1] >=0) # Output power of batteries [MW_el]

        @variable(model, Ph2_electrolyzer[e=1:size(Electrolyzer,1), t=1:timestep] >=0) # Electricicty consumed by each electrolyzer [MW_el]
        @variable(model, P_h2top[h=1:size(H2toP,1), t=1:timestep] >=0) # Electricicty produced by each H2-to-Power technology [MW_el]
        @variable(model, Qh2_dem[n=1:size(Buses,1), t=1:timestep] >=0) # H2 production dedicated for annual demand at each node [MWh_h2/h]
        @variable(model, F_h2[k_h2=1:size(H2PipeLines, 1), t=1:timestep]) # Power flowing through each h2 pipeline [MW]
        @variable(model, Qh2_import[i=1:size(H2Imports,1), t=1:timestep] >=0) # Electricicty produced by H2-to-Power technologies at each node [MW_el]

        @variable(model, Soc_h2[j=1:size(H2Sto,1), t=1:timestep+1] >=0) # State-of-charge of H2 storage at each node [MWh_h2]
        @variable(model, Qh2_sto_in[j=1:size(H2Sto,1), t=1:timestep+1] >=0) # Input of H2 storage at each node [MWh_h2/h]
        @variable(model, Qh2_sto_out[j=1:size(H2Sto,1), t=1:timestep+1] >=0) # Output of H2 storage at each node [MWh_h2/h]

        @variable(model, Sremove[n=1:size(Buses,1), t=1:timestep] >=0) # Energy shifted (removed) at each node [MWh_el/h]
        @variable(model, Sadd[n=1:size(Buses,1), t=1:timestep] >=0) # Energy shifted (added) at each node [MWh_el/h]
        @variable(model, Eshed[n=1:size(Buses,1), t=1:timestep] >=0) #Energy shed at each node [MWh_el/h]

        @variable(model, Pflex_heat[n=1:size(Buses,1), t=1:timestep] >=0) # Heat Flexibility HP1M [MWh_el/h]
        @variable(model, Pflex_mob[n=1:size(Buses,1), t=1:timestep] >=0) # Mobility Flexibility HP1M [MWh_el/h]
        @variable(model, Qflex[n=1:size(Buses,1), t=1:timestep] >=0) # H2 Demand Flexibility [MWh_h2/h]
        #@variable(model, Lf_nu_hot[g=1:size(idx_nu,1)] >= 0)

        ## Constraints 

        # Power Limit of each Generator
        @constraint(model, maxpowerconv[g=1:size(Gen, 1), t=1:timestep], Pconv[g, t] <= availconv[g,t]*pconv_max[g])
        @constraint(model, minpowerconv[g=1:size(Gen, 1), t=1:timestep], Pconv[g, t] >= availconv[g,t]*pconv_min[g])

        # Seasonal nuclear power limits
        #@constraint(model, nukemax[g=1:size(idx_nu, 1), t=1:timestep], Pconv[idx_nu[g], t] == availconv[idx_nu[g],t]*pconv_max[idx_nu[g]])
        #@constraint(model, nukecold1[g=1:size(idx_nu, 1), t=1:2159], Pconv[idx_nu[g], t] == availconv[idx_nu[g],t]*1*pconv_max[idx_nu[g]])
        #@constraint(model, nukehot[g=1:size(idx_nu, 1), t=2160:6551], Pconv[idx_nu[g], t] == availconv[idx_nu[g],t]*Lf_nu_hot[g]*pconv_max[idx_nu[g]])
        #@constraint(model, nukecold2[g=1:size(idx_nu, 1), t=6552:8760], Pconv[idx_nu[g], t] == availconv[idx_nu[g],t]*1*pconv_max[idx_nu[g]])

        # Power Limit of Imports
        @constraint(model, maxpowerimp[i=1:size(Imp, 1), t=1:timestep], Pimp[i, t] <= availimp[i,t]*pimp_max[i] )
        @constraint(model, [t=1:timestep], sum(Pimp[:, t]) <= pimp_max_tot)
        @constraint(model,  [i=1:size(Imp,1)], sum(Pimp[i,t] for t=1:timestep) <= eimp_year_max[i])

        # Power flow Limit through the Lines 
        @constraint(model, maxline1[k=1:size(Lines, 1), t=1:timestep], F[k, t] <= f_max[k]*availline[k,t])
        @constraint(model, maxline2[k=1:size(Lines, 1), t=1:timestep], F[k, t] >= -f_max[k]*availline[k,t])
        
        # Adjacent phase node angle limit
        @constraint(model, [k=1:size(Lines, 1), t=1:timestep], (Theta[l_fromnidx[k], t]-Theta[l_tonidx[k], t]) >= angmin_rad)
        @constraint(model, [k=1:size(Lines, 1), t=1:timestep], (Theta[l_fromnidx[k], t]-Theta[l_tonidx[k], t]) <= angmax_rad)

        # Phase Angle relation
        @constraint(model, pf1[k=1:size(Lines, 1), t=1:timestep], F[k, t] - sb*b_pu[k]*(Theta[l_fromnidx[k], t]-Theta[l_tonidx[k], t])<= (1-availline[k,t])*M)
        @constraint(model, pf2[k=1:size(Lines, 1), t=1:timestep], F[k, t] - sb*b_pu[k]*(Theta[l_fromnidx[k], t]-Theta[l_tonidx[k], t]) >= -(1-availline[k,t])*M)

        #Electrical Power Balance
        @constraint(model, elecbalance[n=1:size(Buses,1), t=1:timestep], sum(Pconv[gn[n], t]) + sum(P_h2top[h2pn[n], t]) + sum(Pimp[impn[n], t]) + sum(Pphs_sto_out[phsn[n], t]) + sum(Pbatt_sto_out[battn[n], t]) + sum(F[l_ton[n], t]) - sum(F[l_fromn[n], t]) - sum(Pphs_sto_in[phsn[n], t]) - sum(Pbatt_sto_in[battn[n], t]) - Curtail[n,t] + Ens[n, t] + Eshed[n, t] - (1+compr_H2[1])*sum(Ph2_electrolyzer[eln[n], t]) - Pflex_heat[n,t] - Pflex_mob[n,t] + Sremove[n,t] - Sadd[n,t] == d_l[n]*load[t] + d_ev[n]*load_ev_fixed[t] + d_hp[n]*load_hp_fixed[t] - d_distgen[n]*pdistgen[t] - d_pv[n]*solar_lf[t]*psolar - d_onshore[n]*onshore_lf[t]*ponshore - d_hydro[n]*hydro_lf[t]*phydro - sum(res[resn[n], t]))
        @constraint(model, limitens[n=1:size(Buses,1), t=1:timestep], Ens[n,t] - Pflex_mob[n,t] - Pflex_heat[n,t] <= d_l[n]*load[t] + d_ev[n]*load_ev_fixed[t] + d_hp[n]*load_hp_fixed[t] )
        @constraint(model, limitcurtail[n=1:size(Buses,1), t=1:timestep], Curtail[n,t] - sum(Pconv[gn[n], t]) <=  d_distgen[n]*pdistgen[t] + d_pv[n]*solar_lf[t]*psolar + d_onshore[n]*onshore_lf[t]*ponshore + d_hydro[n]*hydro_lf[t]*phydro + sum(res[resn[n], t]) )

        #Hub Node, null phase
        @constraint(model, [t=1:timestep], Theta[39, t]==0) 

        #PHS Storage Management
        @constraint(model, [p=1:size(Phs, 1)], Soc_phs[p, 1] ==  eff_pump[p]*Pphs_sto_in[p, 1]  )
        @constraint(model, [p=1:size(Phs, 1)], Pphs_sto_out[p, 1] == 0)
        @constraint(model, [p=1:size(Phs, 1),t=1:timestep], Soc_phs[p, t+1] == Soc_phs[p,t] + eff_pump[p]*Pphs_sto_in[p, t+1]  - Pphs_sto_out[p, t+1]/eff_turb[p])
        @constraint(model, [p=1:size(Phs, 1),t=1:timestep], Pphs_sto_in[p, t] <= p_phs_max[p])
        @constraint(model, [p=1:size(Phs, 1),t=1:timestep], Pphs_sto_out[p, t] <= p_phs_max[p])
        @constraint(model, [p=1:size(Phs, 1),t=1:timestep+1], Soc_phs[p, t] <= soc_phs_max[p])

        #Batteries Management
        @constraint(model, [b=1:size(Batt, 1)], Soc_batt[b, 1] ==  eff_batt_in[b]*Pbatt_sto_in[b, 1]  )
        @constraint(model, [b=1:size(Batt, 1)], Pbatt_sto_out[b, 1] == 0)
        @constraint(model, [b=1:size(Batt, 1),t=1:timestep], Soc_batt[b, t+1] == Soc_batt[b,t] + eff_batt_in[b]*Pbatt_sto_in[b, t+1]  - Pbatt_sto_out[b, t+1]/eff_batt_out[b])
        @constraint(model, [b=1:size(Batt, 1),t=1:timestep], Pbatt_sto_in[b, t] <= p_batt_max[b])
        @constraint(model, [b=1:size(Batt, 1),t=1:timestep], Pbatt_sto_out[b, t] <= p_batt_max[b])
        @constraint(model, [b=1:size(Batt, 1),t=1:timestep+1], Soc_batt[b, t] <= soc_batt_max[b])

        #H2 Energy Balance
        @constraint(model, h2balance[n=1:size(Buses,1), t=1:timestep], Ens_h2[n,t] + sum(Ph2_electrolyzer[eln[n], t])*eff_electrolyzer + sum(Qh2_sto_out[h2ston[n], t]) -  sum(Qh2_sto_in[h2ston[n], t]) + sum(F_h2[h2pipe_ton[n], t]) + sum(Qh2_import[h2impn[n], t]) == Qflex[n,t] + d_h2[n]*load_h2_fixed[t] + sum(F_h2[h2pipe_fromn[n], t]) + sum(P_h2top[h2pn[n], t])/eff_h2top)
        
        # H2 flow Limit through the H2 PipeLines
        @constraint(model, [k_h2=1:size(H2PipeLines, 1), t=1:timestep], F_h2[k_h2, t] <= f_h2_max[k_h2])
        @constraint(model, [k_h2=1:size(H2PipeLines, 1), t=1:timestep], F_h2[k_h2, t] >= -f_h2_max[k_h2])

        # Electrolyzers
        @constraint(model, maxel[e=1:size(Electrolyzer, 1), t=1:timestep], Ph2_electrolyzer[e, t]/ph2_electrolyzer_max[e] <= availel[e,t])
        availel_tot = sum(availel)/timestep/size(Electrolyzer,1)
        println("availel_tot", availel_tot)
        println("min_lf_el", min_lf_el)
        println(sum(ph2_electrolyzer_max)*timestep)
        @constraint(model, minellfavail, sum(Ph2_electrolyzer)/sum(ph2_electrolyzer_max)/timestep >= min_lf_el*availel_tot)

        # Hydrogen-to-power units
        @constraint(model, maxh2p[h=1:size(H2toP, 1), t=1:timestep], P_h2top[h, t] / p_h2top_max[h] <= availh2p[h,t])
        
        # H2 imports
        @constraint(model, [i=1:size(H2Imports,1), t=1:timestep], Qh2_import[i,t]<= qh2_import_max[i])

        #H2 Storage Management
        @constraint(model, [j=1:size(H2Sto,1)], Soc_h2[j,1] ==  Qh2_sto_in[j,1]  )
        @constraint(model, [j=1:size(H2Sto,1)], Qh2_sto_out[j,1] == 0)
        @constraint(model, [j=1:size(H2Sto,1), t=1:timestep], Soc_h2[j,t+1] == Soc_h2[j,t] + eff_h2sto_in[j]*Qh2_sto_in[j,t+1]  - Qh2_sto_out[j,t+1]/eff_h2sto_out[j])
        @constraint(model, [j=1:size(H2Sto,1), t=1:timestep+1], Soc_h2[j,t] <= soc_h2_max[j])
        @constraint(model, [j=1:size(H2Sto, 1),t=1:timestep], Qh2_sto_in[j, t] <= q_h2sto_in_max[j])
        @constraint(model, [j=1:size(H2Sto, 1),t=1:timestep], Qh2_sto_out[j, t] <= q_h2sto_out_max[j])

        # Flexibilities

        temp_heatfex = Vector{JuMP.ConstraintRef}

        #daily h2 and elec flexibility         
        for j=1:Int(ceil(timestep/24))
            start = (j-1)*24+1
            if j>timestep/24
                stop = timestep
            else
                stop = j*24
            end

            # # Daily electrical flexibility from load shifting
            @constraint(model, [n=1:size(Buses,1)], sum(Sremove[n,t] for t=start:stop) <= d_s[n]*sremove_day_max)
            @constraint(model, [n=1:size(Buses,1)], sum(Sremove[n,t] for t=start:stop) == sum(Sadd[n,t] for t=start:stop))
            
            # # Daily electrical flexibility from load sheding
            @constraint(model, [n=1:size(Buses,1)], sum(Eshed[n,t] for t=start:stop) <= d_sh[n]*eshed_day_max)
        
            # # Daily electrical flexibility from heating and mobility sector
            @constraint(model, [n=1:size(Buses,1)], sum(Qflex[n,t] for t=start:stop) == d_h2[n]*load_h2_flex/365)
            @constraint(model, [n=1:size(Buses,1)], sum(Pflex_heat[n,t] for t=start:stop) == d_hp[n]*load_hp_flex[j])
            @constraint(model, [n=1:size(Buses,1)], sum(Pflex_mob[n,t] for t=start:stop) == d_ev[n]*load_ev_flex_daily)

        end

        # Maximum shed power
        @constraint(model, [n=1:size(Buses,1), t=1:timestep], Eshed[n,t] <= d_sh[n]*pshed_max)

        # Maximum shifted power
        @constraint(model, [n=1:size(Buses,1), t=1:timestep], Sremove[n,t] <= d_s[n]*premove_max)

        # Flexible EVs
        @constraint(model, [n=1:size(Buses,1), t=1:timestep], Pflex_mob[n,t] <= d_ev[n]*n_ev*w_v1m*avail_v1m[t]*p_ev_charger)

        # Flexible HPs
        @constraint(model, limitheatflex[n=1:size(Buses,1), t=1:timestep], Pflex_heat[n,t] <= d_hp[n]*n_hp*w_hp1m*p_hp_individual)

        # Objective function 
        if ((size(Imp)[1] > 0) && (size(H2Imports)[1] > 0)) # If interconnections
            println("imp + h2imp")
            @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(c_imp.*Pimp) + sum(ch2imp.*Qh2_import))
        elseif ((size(Imp)[1] == 0) && (size(H2Imports)[1] > 0))
            println("h2imp")
            @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(ch2imp.*Qh2_import))
        elseif ((size(Imp)[1] > 0) && (size(H2Imports)[1] == 0))
            println("imp")
            @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(c_imp.*Pimp))
        else
            println("non")
            @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2) + cshift*sum(Sremove) + cshed*sum(Eshed) )
        end

        stop_modelopt = now()
        runtime_modelopt = Dates.value(stop_modelopt - start_modelopt)/1000 #milliseconds to seconds
        println("--- Thread ", Threads.threadid(),"- Time : ",hour(stop_modelopt), ":", minute(stop_modelopt), ":", second(stop_modelopt), " - End Optimization model : ", runtime_modelopt, " seconds") # ~13.8 s

## Batch loop ##

        for batch in 1:n_batch

            # Threads.lock(u) do
            #     global n_mc = n_mc + n_years_in_batch
            #     println("\n", "thread ", Threads.threadid(), " - n_mc= ", n_mc) # ~13.8 s
            # end

## Initialization of batch variables

            lol_set_batch = Vector{Int64}(undef, n_years_in_batch)
            lol_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)
            loe_set_batch = Vector{Float64}(undef, n_years_in_batch)
            loe_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)
            h2loe_set_batch = Vector{Float64}(undef, n_years_in_batch)
            h2loe_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)

            e_curt_set_batch = Vector{Float64}(undef, n_years_in_batch)
            e_curt_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)
            p_curt_set_batch = Vector{Float64}(undef, n_years_in_batch)

            phscycle_set_batch = Array{Float64}(undef, size(Phs)[1], n_years_in_batch)
            battcycle_set_batch = Array{Float64}(undef, size(Batt)[1], n_years_in_batch)
            h2cycle_set_batch = Array{Float64}(undef, size(H2Sto)[1], n_years_in_batch)

            shift_set_batch = Vector{Float64}(undef, n_years_in_batch)
            shift_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)
            shed_set_batch = Vector{Float64}(undef, n_years_in_batch)
            shed_nodal_set_batch = Array{Float64}(undef, size(Buses)[1], n_years_in_batch)

            conv_set_batch = Vector{Float64}(undef, n_years_in_batch)
            genlf_set_batch = Array{Float64}(undef, size(Gen)[1], n_years_in_batch)
            offshore_set_batch = Vector{Float64}(undef, n_years_in_batch)

            imp_set_batch = Vector{Float64}(undef, n_years_in_batch)
            imp_nodal_set_batch = Array{Float64}(undef, size(Imp)[1], n_years_in_batch)
            h2imp_set_batch = Vector{Float64}(undef, n_years_in_batch) 
            h2imp_nodal_set_batch = Array{Float64}(undef, size(H2Imports)[1], n_years_in_batch)

            lly_set_batch = Array{Float64}(undef, size(Lines)[1], n_years_in_batch)
            plly_set_batch = Array{Float64}(undef, size(H2PipeLines)[1], n_years_in_batch)

            elinst_set_batch = Vector{Float64}(undef, n_years_in_batch)
            elinst_nodal_set_batch = Array{Float64}(undef, size(Electrolyzer)[1], n_years_in_batch)
            ellf_set_batch = Vector{Float64}(undef, n_years_in_batch) 
            ellf_nodal_set_batch = Array{Float64}(undef, size(Electrolyzer)[1], n_years_in_batch)

            h2topinst_set_batch = Vector{Float64}(undef, n_years_in_batch)  
            h2topinst_nodal_set_batch = Array{Float64}(undef, size(H2toP)[1], n_years_in_batch)
            h2toplf_set_batch = Vector{Float64}(undef, n_years_in_batch) 
            h2toplf_nodal_set_batch = Array{Float64}(undef, size(H2toP)[1], n_years_in_batch)

            soc_h2_max_res_set_batch = Vector{Float64}(undef, n_years_in_batch) 

            cost_set_batch = Vector{Float64}(undef, n_years_in_batch)
            runtimeopti_set_batch = Vector{Float64}(undef, n_years_in_batch)
## Monte-Carlo run

            for year in 1:n_years_in_batch
                println("--- Thread ", Threads.threadid()," - Time : ", hour(now()), ":", minute(now()), ":", second(now()),  " - New Monte Carlo, year ", year, " in batch ", batch)

## 4. Update uncertain parameters ##

                start_updateparam = now()
                println("--- Thread ", Threads.threadid()," - Time : ", hour(start_updateparam), ":", minute(start_updateparam), ":", second(start_updateparam),  " - Start updating param")

                # Availabilities due to potential outages
                if outagestech == true

                    availconv_update = avail_comp(length(pconv_max), timestep, mttf_conv, mttr_conv, pconv_max, false, zeros(timestep), seed1)[1] # Availability of generators (TTF and TTR sampling)
                    availline_update = avail_comp(lastindex(f_max), timestep, mttf_f, mttr_f, f_max, false, zeros(timestep), seed2)[1] # Availability of Lines (TTF and TTR sampling)

                    if size(Electrolyzer)[1] > 0
                        availel_update = avail_comp(lastindex(ph2_electrolyzer_max), timestep, mttf_el, mttr_el, ph2_electrolyzer_max, false, zeros(timestep), seed3)[1] # Availability of electrolyzers (TTF and TTR sampling)
                    end

                    if size(H2toP)[1] > 0
                    availh2p_update = avail_comp(lastindex(p_h2top_max), timestep, mttf_h2p, mttr_h2p, p_h2top_max, false, zeros(timestep), seed4)[1] # Availability of hydrogen-to-power units (TTF and TTR sampling)
                    end

                else
                    availconv_update = ones(length(pconv_max), timestep)
                    availline_update = ones(lastindex(f_max), timestep)
                    availel_update = ones(lastindex(ph2_electrolyzer_max), timestep)
                    availh2p_update = ones(lastindex(p_h2top_max), timestep)
                    println("availel_tot_update", availel_tot_update)
                end

                # Availabilities of interconnections and neighbouring generation potential
                if size(Imp)[1] > 0 # If interconnections
                    if import_avail == "random"
                        availimp_update = 2*rand(size(Imp,1), timestep)
                        replace(x -> x>1 ? 0 : x, availimp_update)
                        c_imp_update = c_imp
                    elseif import_avail == "max"
                        availimp_update = avail_comp(lastindex(pimp_max), timestep, mttf_imp, mttr_imp, pimp_max, false, zeros(timestep), seed5)[1]
                        c_imp_update = Array{Float64}(undef, size(import_countries)[1], timestep)
                        for c in 1:size(import_countries)[1]
                            Random.seed!(seed_noise_c_imp[c]+n_mc) # Setting the seed
                            dist_c_imp = Normal(mean_noise_c_imp[c], std_noise_c_imp[c])
                            noise = rand(dist_c_imp, (timestep))
                            c_imp_update[c, :] = c_imp_base[c, :] + noise
                        end
                    end   
                else
                    c_imp_update = 0
                    availimp_update = zeros(size(Imp,1), timestep)
                end

                # Offshore wind
                if size(Res)[1] > 0 
                    res_update_f = offshore_time_series4(nb_wt_offshore, wf_names, timestep, mttf_res, mttr_res, outageswind, constance, seed6)
                    res_update = increasefactoroffshore*res_update_f[1]
                    r_wind_update = res_update_f[2]
                else
                    res_update = zeros(2,timestep)
                end

                #Load
                load_update = load_elia_pattern*yearly_demand_elec

                availel_tot_update = sum(availel_update)/timestep/size(Electrolyzer,1)
                println("availel_tot", availel_tot_update)
                stop_updateparam = now()
                runtime_updateparam = Dates.value(stop_updateparam - start_updateparam)/1000 #milliseconds to seconds
                println("--- Thread ", Threads.threadid(), " - Time : ", hour(stop_updateparam), ":", minute(stop_updateparam), ":", second(stop_updateparam), " - End update parameters: ", runtime_updateparam, " seconds")


## 5. Update constraints ##
                start_updateconstr = now()
                println("--- Thread ", Threads.threadid()," - Time : " ,hour(stop_updateparam), ":", minute(stop_updateparam), ":", second(stop_updateparam),  " - Start update constraints")

                for t=1:timestep
                    for g in eachindex(pconv_max)
                        set_normalized_rhs(maxpowerconv[g,t], availconv_update[g,t]*pconv_max[g])
                        set_normalized_rhs(minpowerconv[g,t], availconv_update[g,t]*pconv_min[g])
                    end

                    for i in eachindex(pimp_max)
                        set_normalized_rhs(maxpowerimp[i,t], availimp_update[i,t]*pimp_max[i])
                    end

                    for n=1:size(Buses,1)
                        set_normalized_rhs(elecbalance[n,t], d_l[n]*load_update[t] + d_ev[n]*load_ev_fixed[t] + d_hp[n]*load_hp_fixed[t] - d_distgen[n]*pdistgen[t] - d_pv[n]*solar_lf[t]*psolar - d_onshore[n]*onshore_lf[t]*ponshore - d_hydro[n]*hydro_lf[t]*phydro - sum(res_update[resn[n], t]))
                        set_normalized_rhs(limitens[n,t], d_l[n]*load_update[t] + d_ev[n]*load_ev_fixed[t] + d_hp[n]*load_hp_fixed[t])
                        set_normalized_rhs(limitcurtail[n,t], d_distgen[n]*pdistgen[t] + d_pv[n]*solar_lf[t]*psolar + d_onshore[n]*onshore_lf[t]*ponshore + d_hydro[n]*hydro_lf[t]*phydro + sum(res_update[resn[n], t]) )
                    end

                    for k in eachindex(f_max)
                        set_normalized_rhs(maxline1[k,t], availline_update[k,t]*f_max[k])
                        set_normalized_rhs(maxline2[k,t], -availline_update[k,t]*f_max[k])

                        set_normalized_rhs(pf1[k,t], (1-availline_update[k,t])*M)
                        set_normalized_rhs(pf2[k,t], -(1-availline_update[k,t])*M)
                    end

                    for h in eachindex(ph2_electrolyzer_max)
                        set_normalized_rhs(maxel[h,t], availel_update[h,t])
                    end

                    for h in eachindex(p_h2top_max)
                        set_normalized_rhs(maxh2p[h,t], availh2p_update[h,t])
                    end
                end

                set_normalized_rhs(minellfavail, min_lf_el*availel_tot_update)


                if ((size(Imp)[1] > 0) && (size(H2Imports)[1] > 0)) # If interconnections
                    @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(c_imp_update.*Pimp) + sum(ch2imp.*Qh2_import))
                elseif ((size(Imp)[1] == 0) && (size(H2Imports)[1] > 0))
                    @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(ch2imp.*Qh2_import))
                elseif ((size(Imp)[1] > 0) && (size(H2Imports)[1] == 0))
                    @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2)  + cshift*sum(Sremove) + cshed*sum(Eshed) + sum(c_imp_update.*Pimp))
                else
                    @objective(model, Min, sum(c_gen_op.*Pconv) + voll*sum(Ens) + voll_h2*sum(Ens_h2) + cshift*sum(Sremove) + cshed*sum(Eshed) )
                end


                stop_updateconstr = now()
                runtime_updateconstr = Dates.value(stop_updateconstr - start_updateconstr)/1000 #milliseconds to seconds
                println("--- Thread ", Threads.threadid()," - Time : " ,hour(stop_updateconstr), ":", minute(stop_updateconstr), ":", second(stop_updateconstr)," - End update constraints: ", runtime_updateconstr, " seconds")

## 6. Optimize updated problem##
                MathOptInterface.Utilities.reset_optimizer(model) #Reset optimizer to start searching for new solution NOT from previous solution
                start_opti = now()
                println("--- Thread ", Threads.threadid(), " - Time : ",hour(start_opti), ":", minute(start_opti), ":", second(start_opti), " - Start optimization")
                optimize!(model)
                stop_opti = now()
                runtime_opti = Dates.value(stop_opti - start_opti)/1000/60 #milliseconds to minutes
                println("--- Thread ", Threads.threadid()," - Time : ",hour(stop_opti), ":", minute(stop_opti), ":", second(stop_opti), " - End optimization: ", round(runtime_opti, digits=4), " minutes", )

## 7. Print Inputs
                ## Print input values: if optimization failed OR high lol OR first years

                if isempty(lole)
                    #println("The lole is empty")
                    last_lole = 0
                else
                    last_lole = last(lole)
                end

                start_writecompleteyear = now()

                println("--- Thread ", Threads.threadid()," - Time : " ,hour(start_writecompleteyear), ":", minute(start_writecompleteyear), ":", second(start_writecompleteyear), " - Start writing")
                
                #Pint input dynamic paramters
                #If model did not converge - if LOL higher than 10*LOLE - if first year of batch 1
                if has_values(model) == false || (length(get_index_nonzeros_multinode(value.(Ens))) >= 10*(1+last_lole) && n_mc > 60) || ((year == 1) && batch == 1 && print_yearlydata == 1)
                    if has_values(model) == false 
                        #println("Optimization failed")
                        printfolder = "failuremc"
                    elseif (length(get_index_nonzeros_multinode(value.(Ens))) >= 10*(1+last_lole) && n_mc > 60)
                        #println("Unexpected high LOL")
                        printfolder = "highLOL"
                    else
                        #println("Print first years")
                        printfolder = "firstyears"
                    end

                    println("--- Thread ", Threads.threadid()," - Time : " ,hour(start_writecompleteyear), ":", minute(start_writecompleteyear), ":", second(start_writecompleteyear)," - Start writing in ",printfolder )
                    

                    #Availabilities of generators
                    open("./Results/"*simu_name*"/"*printfolder*"/avail_conv"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            writedlm(io, [availconv_update[:,t]])
                        end
                    end

                    #Availabilities of lines
                    open("./Results/"*simu_name*"/"*printfolder*"/avail_lines"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            writedlm(io, [availline_update[:,t]])
                        end
                    end

                    #Availabilities of imports
                    open("./Results/"*simu_name*"/"*printfolder*"/avail_imp"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            writedlm(io, [availimp_update[:,t]])
                        end
                    end

                    # #Availabilities of electrolyzers
                    # open("./Results/"*simu_name*"/"*printfolder*"/avail_el"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                    #     for t=1:timestep
                    #         writedlm(io, [availel_update[:,t]])
                    #     end
                    # end

                    # #Availabilities of h2p units
                    # open("./Results/"*simu_name*"/"*printfolder*"/avail_h2p"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                    #     for t=1:timestep
                    #         writedlm(io, [availh2p_update[:,t]])
                    #     end
                    # end

                    # #Load
                    # open("./Results/"*simu_name*"/"*printfolder*"/load"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                    #     for t=1:timestep
                    #         writedlm(io, [load_update[t]])
                    #     end
                    # end

                    #Res
                    open("./Results/"*simu_name*"/"*printfolder*"/res"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        # for t=1:timestep
                        #     writedlm(io, [res_update[:,t]])
                        # end
                        writedlm(io,r_wind_update)
                    end
                end

## 8. Print outputs 

                #Pint dtailed output data
                #If LOL higher than 10*LOLE - if first year of batch 1
                if (year == 1 && batch == 1 && print_yearlydata == 1) || (length(get_index_nonzeros_multinode(value.(Ens))) >= 10*(1+last_lole) && n_mc > 60)
                    if (year == 1) && batch == 1 && print_yearlydata == 1
                        printfolder = "firstyears"
                    else
                        printfolder = "highLOL"
                    end
                    
                    # Elec Balance: [Buses, Timestep]
                    open("./Results/"*simu_name*"/"*printfolder*"/elecbalance"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for n=1:size(Buses,1)
                            for t=1:timestep
                                if length(gn[n])==0
                                    output_pconv = 0
                                else
                                    output_pconv = sum(value.(Pconv[gn[n], t])[x,:] for x=1:length(gn[n]))
                                end

                                if length(resn[n])==0
                                    output_res = 0
                                else
                                    output_res = sum(value.(res_update[resn[n], t])[x,:] for x=1:length(resn[n]))
                                end

                                if length(impn[n])==0
                                    output_imp = 0
                                else
                                    output_imp = sum(value.(Pimp[impn[n], t])[x,:] for x=1:length(impn[n]))
                                end

                                if length(phsn[n])==0
                                    output_phs_out = 0
                                    output_phs_in = 0
                                else
                                    output_phs_out = sum(value.(Pphs_sto_out[phsn[n], t])[x,:] for x=1:length(phsn[n]))
                                    output_phs_in = sum(value.(Pphs_sto_in[phsn[n], t])[x,:] for x=1:length(phsn[n]))
                                end

                                if length(battn[n])==0
                                    output_batt_out = 0
                                    output_batt_in = 0
                                else
                                    output_batt_out = sum(value.(Pbatt_sto_out[battn[n], t])[x,:] for x=1:length(battn[n]))
                                    output_batt_in = sum(value.(Pbatt_sto_in[battn[n], t])[x,:] for x=1:length(battn[n]))
                                end

                                if length(l_ton[n])==0
                                    output_fto = 0
                                else
                                    output_fto = sum(value.(F[l_ton[n], t])[x,:] for x=1:length(l_ton[n]))
                                end

                                if length(l_fromn[n])==0
                                    output_ffrom = 0
                                else
                                    output_ffrom = sum(value.(F[l_fromn[n], t])[x,:] for x=1:length(l_fromn[n]))
                                end

                                if length(eln[n])==0
                                    output_el = 0
                                else
                                    output_el = sum(value.(Ph2_electrolyzer[eln[n], t])[x,:] for x=1:length(eln[n]))
                                end

                                if length(h2pn[n])==0
                                    output_h2top = 0
                                else
                                    output_h2top = sum(value.(P_h2top[h2pn[n], t])[x,:] for x=1:length(h2pn[n]))
                                end
                                    
                                output_ens = value(Ens[n,t])
                                output_load = d_l[n]*load_update[t]
                                output_curtail = value(Curtail[n,t])
                                output_eshed = value(Eshed[n,t])
                                output_ev_fixed = d_ev[n]*load_ev_fixed[t]
                                output_hp_fixed = d_hp[n]*load_hp_fixed[t]
                                output_ev_flex = value(Pflex_mob[n,t])
                                output_hp_flex = value(Pflex_heat[n,t])
                                output_sremove = value(Sremove[n,t])
                                output_sadd = value(Sadd[n,t])
                                output_pv = d_pv[n]*solar_lf[t]*psolar
                                output_onshore = d_onshore[n]*onshore_lf[t]*ponshore
                                output_hydro = d_hydro[n]*hydro_lf[t]*phydro
                                output_distgen = d_distgen[n]*pdistgen[t]


                                writedlm(io, [n t output_pconv output_distgen output_res output_fto output_imp output_phs_out output_batt_out output_ens output_load output_el compr_H2[1]*output_el output_h2top output_ffrom output_phs_in output_batt_in output_curtail output_eshed output_ev_fixed output_ev_flex output_hp_fixed output_hp_flex output_sremove output_sadd output_pv output_onshore output_hydro])
                            end
                        end      
                    end
                    
                    # H2 Balance
                    open("./Results/"*simu_name*"/"*printfolder*"/h2balance"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        
                        for n=1:size(Buses,1)
                            for t=1:timestep

                                if length(h2pipe_ton[n])==0
                                    output_fh2to = 0
                                else
                                    output_fh2to = sum(value.(F_h2[h2pipe_ton[n], t])[x,:] for x=1:length(h2pipe_ton[n]))
                                end

                                if length(h2pipe_fromn[n])==0
                                    output_fh2from = 0
                                else
                                    output_fh2from = sum(value.(F_h2[h2pipe_fromn[n], t])[x,:] for x=1:length(h2pipe_fromn[n]))
                                end

                                if length(eln[n])==0
                                    output_el = 0
                                else
                                    output_el = sum(value.(Ph2_electrolyzer[eln[n], t])[x,:] for x=1:length(eln[n]))
                                end

                                if length(h2pn[n])==0
                                    output_h2top = 0
                                else
                                    output_h2top = sum(value.(P_h2top[h2pn[n], t])[x,:] for x=1:length(h2pn[n]))
                                end

                                if length(h2ston[n])==0
                                    output_h2_sto_in = 0
                                    output_h2_sto_out = 0
                                else
                                    output_h2_sto_in = sum(value.(Qh2_sto_in[h2ston[n], t])[x,:] for x=1:length(h2ston[n]))
                                    output_h2_sto_out = sum(value.(Qh2_sto_out[h2ston[n], t])[x,:] for x=1:length(h2ston[n]))
                                end

                                if length(h2impn[n])==0
                                    output_h2import = 0
                                else
                                    output_h2import = sum(value.(Qh2_import[h2impn[n], t])[x,:] for x=1:length(h2impn[n]))
                                end


                                output_ensh2 = value(Ens_h2[n,t])
                                output_qh2fixed = d_h2[n]*load_h2_fixed[t]
                                output_qh2flex = value(Qflex[n,t])

                                writedlm(io, [n t output_ensh2 output_h2import output_el output_el*eff_electrolyzer output_fh2to output_fh2from output_qh2fixed output_qh2flex output_h2top/eff_h2top output_h2_sto_out output_h2_sto_in])
                            end
                        end
                    end

                    #Line Power Flows
                    open("./Results/"*simu_name*"/"*printfolder*"/lines_power"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                                output_alline_power_t = value.(F[:,t])
                                writedlm(io, [output_alline_power_t])
                        end
                    end

                    #Line Loadings 
                    open("./Results/"*simu_name*"/"*printfolder*"/lines_loading"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            output_alline_loading_t = value.(F[:,t])./f_max
                            writedlm(io, [output_alline_loading_t])
                        end
                    end

                    #Load factors of generators
                    open("./Results/"*simu_name*"/"*printfolder*"/genlf"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            writedlm(io, [value.(Pconv[:,t])./pconv_max])
                        end
                    end
                end

                stop_writecompleteyear = now()
                runtime_writecompleteyear = Dates.value(stop_writecompleteyear - start_writecompleteyear)/1000/60 #milliseconds to minutes
                println("--- Thread ", Threads.threadid()," - Time : ",hour(stop_writecompleteyear), ":", minute(stop_writecompleteyear), ":", second(stop_writecompleteyear), " - Stop writng: ", round(runtime_writecompleteyear, digits=4), " minutes")


                
                ## Verify elecbalance if first year in batch 1
                if (year == 1) && batch == 1 && print_verifyelec == 1
                    start_verifelecbalance = now()
                    println("--- Thread ", Threads.threadid()," - Time : ",hour(start_verifelecbalance), ":", minute(start_verifelecbalance), ":", second(start_verifelecbalance), " - Start elec balance verification")
                    for n=1:size(Buses,1)
                        for t=1:timestep
                            eq = abs(value(sum(Pconv[gn[n], t]))+ d_distgen[n]*pdistgen[t] + value(sum(P_h2top[h2pn[n], t])) + value(sum(Pimp[impn[n], t])) + value(sum(Pphs_sto_out[phsn[n], t])) + value(sum(Pbatt_sto_out[battn[n], t])) + value(sum(F[l_ton[n], t])) - value(sum(F[l_fromn[n], t])) - value(sum(Pphs_sto_in[phsn[n], t])) - value(sum(Pbatt_sto_in[battn[n], t])) - value(Curtail[n,t]) + value(Ens[n, t]) + value(Eshed[n, t]) - (1+compr_H2[1])*value(sum(Ph2_electrolyzer[eln[n], t])) - d_ev[n]*load_ev_fixed[t] - d_hp[n]*load_hp_fixed[t] - value(Pflex_heat[n, t]) - value(Pflex_mob[n, t]) + value(Sremove[n,t]) - value(Sadd[n,t]) + d_pv[n]*solar_lf[t]*psolar + d_onshore[n]*onshore_lf[t]*ponshore + d_hydro[n]*hydro_lf[t]*phydro - d_l[n]*load_update[t]+sum(res_update[resn[n], t])) 
			                if  eq <= 1e-10
                                #println(eq)
                            else
                                print("p", eq)
                            end 
                        end
                    end
                    stop_verifelecbalance = now()
                    runtime_verifelecbalance = Dates.value(stop_verifelecbalance - start_verifelecbalance)/1000 #milliseconds to seconds
                    println("--- Thread ", Threads.threadid()," - Time : ",hour(stop_verifelecbalance), ":", minute(stop_verifelecbalance), ":", second(stop_verifelecbalance), " - End verification of elec balance : ", runtime_verifelecbalance, " seconds")
                end


                # Print elecbalance, h2 blance, and gen power when ens occurs: 

                if sum(value.(Ens)) > 0 #If ENS occurs: create txt file
                    open("./Results/"*simu_name*"/ens/elecbalance"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for n=1:size(Buses,1)
                            for t=1:timestep

                                if sum(value.(Ens[:,t])) > 0 #If ENS: print 

                                    if length(gn[n])==0
                                        output_pconv = 0
                                    else
                                        output_pconv = sum(value.(Pconv[gn[n], t])[x,:] for x=1:length(gn[n]))
                                    end

                                    if length(resn[n])==0
                                        output_res = 0
                                    else
                                        output_res = sum(value.(res_update[resn[n], t])[x,:] for x=1:length(resn[n]))
                                    end

                                    if length(impn[n])==0
                                        output_imp = 0
                                    else
                                        output_imp = sum(value.(Pimp[impn[n], t])[x,:] for x=1:length(impn[n]))
                                    end

                                    if length(phsn[n])==0
                                        output_phs_out = 0
                                        output_phs_in = 0
                                    else
                                        output_phs_out = sum(value.(Pphs_sto_out[phsn[n], t])[x,:] for x=1:length(phsn[n]))
                                        output_phs_in = sum(value.(Pphs_sto_in[phsn[n], t])[x,:] for x=1:length(phsn[n]))
                                    end

                                    if length(battn[n])==0
                                        output_batt_out = 0
                                        output_batt_in = 0
                                    else
                                        output_batt_out = sum(value.(Pbatt_sto_out[battn[n], t])[x,:] for x=1:length(battn[n]))
                                        output_batt_in = sum(value.(Pbatt_sto_in[battn[n], t])[x,:] for x=1:length(battn[n]))
                                    end

                                    if length(l_ton[n])==0
                                        output_fto = 0
                                    else
                                        output_fto = sum(value.(F[l_ton[n], t])[x,:] for x=1:length(l_ton[n]))
                                    end

                                    if length(l_fromn[n])==0
                                        output_ffrom = 0
                                    else
                                        output_ffrom = sum(value.(F[l_fromn[n], t])[x,:] for x=1:length(l_fromn[n]))
                                    end

                                    if length(eln[n])==0
                                        output_el = 0
                                    else
                                        output_el = sum(value.(Ph2_electrolyzer[eln[n], t])[x,:] for x=1:length(eln[n]))
                                    end

                                    if length(h2pn[n])==0
                                        output_h2top = 0
                                    else
                                        output_h2top = sum(value.(P_h2top[h2pn[n], t])[x,:] for x=1:length(h2pn[n]))
                                    end
                                        
                                    output_ens = value(Ens[n,t])
                                    output_load = d_l[n]*load_update[t]
                                    output_curtail = value(Curtail[n,t])
                                    output_eshed = value(Eshed[n,t])
                                    output_ev_fixed = d_ev[n]*load_ev_fixed[t]
                                    output_hp_fixed = d_hp[n]*load_hp_fixed[t]
                                    output_ev_flex = value(Pflex_mob[n,t])
                                    output_hp_flex = value(Pflex_heat[n,t])
                                    output_sremove = value(Sremove[n,t])
                                    output_sadd = value(Sadd[n,t])
                                    output_pv = d_pv[n]*solar_lf[t]*psolar
                                    output_onshore = d_onshore[n]*onshore_lf[t]*ponshore
                                    output_hydro = d_hydro[n]*hydro_lf[t]*phydro
                                    output_distgen = d_distgen[n]*pdistgen[t]


                                    writedlm(io, [n t output_pconv output_distgen output_res output_fto output_imp output_phs_out output_batt_out output_ens output_load output_el compr_H2[1]*output_el output_h2top output_ffrom output_phs_in output_batt_in output_curtail output_eshed output_ev_fixed output_ev_flex output_hp_fixed output_hp_flex output_sremove output_sadd output_pv output_onshore output_hydro])
                            
                                end
                            end
                        end   
                    end

                    # H2 Balance
                    open("./Results/"*simu_name*"/ens/h2balance"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        
                        for n=1:size(Buses,1)
                            for t=1:timestep
                                if sum(value.(Ens[:,t])) > 0 #If ENS: print 

                                    if length(h2pipe_ton[n])==0
                                        output_fh2to = 0
                                    else
                                        output_fh2to = sum(value.(F_h2[h2pipe_ton[n], t])[x,:] for x=1:length(h2pipe_ton[n]))
                                    end

                                    if length(h2pipe_fromn[n])==0
                                        output_fh2from = 0
                                    else
                                        output_fh2from = sum(value.(F_h2[h2pipe_fromn[n], t])[x,:] for x=1:length(h2pipe_fromn[n]))
                                    end

                                    if length(eln[n])==0
                                        output_el = 0
                                    else
                                        output_el = sum(value.(Ph2_electrolyzer[eln[n], t])[x,:] for x=1:length(eln[n]))
                                    end

                                    if length(h2pn[n])==0
                                        output_h2top = 0
                                    else
                                        output_h2top = sum(value.(P_h2top[h2pn[n], t])[x,:] for x=1:length(h2pn[n]))
                                    end

                                    if length(h2ston[n])==0
                                        output_h2_sto_in = 0
                                        output_h2_sto_out = 0
                                    else
                                        output_h2_sto_in = sum(value.(Qh2_sto_in[h2ston[n], t])[x,:] for x=1:length(h2ston[n]))
                                        output_h2_sto_out = sum(value.(Qh2_sto_out[h2ston[n], t])[x,:] for x=1:length(h2ston[n]))
                                    end

                                    if length(h2impn[n])==0
                                        output_h2import = 0
                                    else
                                        output_h2import = sum(value.(Qh2_import[h2impn[n], t])[x,:] for x=1:length(h2impn[n]))
                                    end


                                    output_ensh2 = value(Ens_h2[n,t])
                                    output_qh2fixed = d_h2[n]*load_h2_fixed[t]
                                    output_qh2flex = value(Qflex[n,t])

                                    writedlm(io, [n t output_ensh2 output_h2import output_el output_el*eff_electrolyzer output_fh2to output_fh2from output_qh2fixed output_qh2flex output_h2top/eff_h2top output_h2_sto_out output_h2_sto_in])
                                end
                            end
                        end
                    end

                    #Load factors of generators
                    open("./Results/"*simu_name*"/ens/genlf"*string(year)*"_batch"*string(batch)*"_thread"*string(Threads.threadid())*".txt", "w") do io
                        for t=1:timestep
                            if sum(value.(Ens[:,t])) > 0
                                writedlm(io, [vcat(t,value.(Pconv[:,t])./pconv_max)])
                            end
                        end
                    end
                end

## Extract values of interest from variables

                lol = length(get_index_nonzeros_multinode(value.(Ens))) # Loss-of-load [hours/year]
                lol_nodal = get_index_nonzeros_2d(value.(Ens)) # Nodal Loss-of-load [hours/year]
                loe = sum(value.(Ens)) # Loss-of-energy [MWh/year]
                loe_nodal = sum(value.(Ens[:,x]) for x=1:timestep) # Nodal Loss-of-energy [MWh/year]
		        h2loe = sum(value.(Ens_h2)) # Loss-of-h2-energy [MWh_h2/year]
                h2loe_nodal = sum(value.(Ens_h2[:,x]) for x=1:timestep) # Nodal Loss-of-h2-energy [MWh_h2/year]

                e_curt = sum(value.(Curtail)) # Energy curtailed [MWh/year]
                e_curt_nodal = sum(value.(Curtail[:,x]) for x=1:timestep) # Nodal Energy curtailed [MWh/year]
                p_curt = maximum(sum(value.(Curtail[x,:]) for x=1:size(Buses,1))) # Peak curtailment [MW]

                phscycle = sum(value.(Pphs_sto_in[:,x]) for x=1:timestep).*eff_pump[:]./soc_phs_max[:] # Number of PHS cycles [-]
                battcycle = sum(value.(Pbatt_sto_in[:,x]) for x=1:timestep).*eff_batt_in[:]./soc_batt_max[:] # Number of battry cycles [-]
                h2cycle = sum(value.(Qh2_sto_in[:,x]) for x=1:timestep)./soc_h2_max[:] # Number of h2 storage units cycles [-]

                shift = sum(value.(Sremove)) # Quantity of nergy shifted [MWh/year]
                shift_nodal = sum(value.(Sremove[:,x]) for x=1:timestep) # Nodal Quantity of nergy shifted [MWh/year]
                shed = sum(value.(Eshed)) # Quantity of energy shed [MWh//year]
                shed_nodal = sum(value.(Eshed[:,x]) for x=1:timestep) # Nodal Quantity of energy shed [MWh//year]

                conv = sum(value.(Pconv)) # Production of conventional units [MWh/year]
                genlf = sum(value.(Pconv[:,x]) for x=1:timestep)./pconv_max/timestep # Individual load factor of conventional units [MWh/year]
                offshore = sum(res_update) # Offshore wind farms production [MWh/year]

                imp = sum(value.(Pimp)) # Total Energy imported [MWh/year]
                imp_nodal = sum(value.(Pimp[:,x]) for x=1:timestep) # Energy imported at each interconnection [MWh/year]
                h2imp = sum(value.(Qh2_import)) # Total H2 Energy imported [MWh/year]
                h2imp_nodal = sum(value.(Qh2_import[:,x]) for x=1:timestep) # H2 Energy imported at each import hub [MWh/year]

                lly = sum(abs.(value.(F[:,x])) for x=1:timestep)./f_max/timestep # Yearly loading of each line [-]
                plly = sum(abs.(value.(F_h2[:,x])) for x=1:timestep)./f_h2_max/timestep # Yearly loading of each pipeline [-]

                if size(Electrolyzer, 1) > 0
                    elinst = maximum(sum(value.(Ph2_electrolyzer[x,:]) for x=1:size(Electrolyzer,1))) # Total Maximum power consumed by electrolyzers [MW]
                    #print("elinst", elinst, size(elinst))
                    ellf = sum(value.(Ph2_electrolyzer))/sum(ph2_electrolyzer_max)/timestep # Total load factor of electrolyzers [-]
                else  
                    elinst = 0 
                    ellf = 0   
                end 
                elinst_nodal = maximum(value.(Ph2_electrolyzer[:,x]) for x=1:timestep) # Individual maximum power consumed by electrolyzers [MW]
                ellf_nodal = sum(value.(Ph2_electrolyzer[:,x]) for x=1:timestep)./ph2_electrolyzer_max/timestep # Individual load factor of electrolyzers [-]

                if size(H2toP, 1) > 0
                    h2topinst = maximum(sum(value.(P_h2top[x,:]) for x=1:size(H2toP,1))) # Total Maximum power produced by hydrogen-to-power units [MW]
                    #println("h2topinst", h2topinst)
                    h2toplf = mean(sum(value.(P_h2top[x,:])/sum(p_h2top_max)  for x=1:size(H2toP,1))) # Total load factor of hydrogen-to-power units [-]
                else    
                    h2topinst = 0 
                    h2toplf = 0
                end
                h2topinst_nodal = maximum(value.(P_h2top[:,x]) for x=1:timestep) # Individual maximum power produced by hydrogen-to-power units [MW]
                h2toplf_nodal = sum(value.(P_h2top[:,x]) for x=1:timestep)./p_h2top_max/timestep # Individual load factors of hydrogen-to-power units [-]
                
                if size(H2Sto, 1) > 0
                    soc_h2_max_res = maximum(sum(value.(Soc_h2[x,:]) for x=1:size(H2Sto,1))) # Total Maximum energy stored in h2 storage units [MWh]
                else
                    soc_h2_max_res = 0
                end
            
                cost = objective_value(model) # Cost [€/year]

                ## Compare some values of interest

## Store values of interest in batch variables

                lol_set_batch[year] = lol
                lol_nodal_set_batch[:, year] = lol_nodal
                loe_set_batch[year] = loe
                loe_nodal_set_batch[:, year] = loe_nodal
                h2loe_set_batch[year] = h2loe
                h2loe_nodal_set_batch[:, year] = h2loe_nodal

                e_curt_set_batch[year] = e_curt
                e_curt_nodal_set_batch[:, year] = e_curt_nodal
                p_curt_set_batch[year] = p_curt

                phscycle_set_batch[:, year] = phscycle
                battcycle_set_batch[:, year] = battcycle
                h2cycle_set_batch[:, year] = h2cycle
                
                shift_set_batch[year] = shift
                shift_nodal_set_batch[:, year] = shift_nodal
                shed_set_batch[year] = shed
                shed_nodal_set_batch[:, year] = shed_nodal

                conv_set_batch[year] = conv
                genlf_set_batch[:,year] = genlf
                offshore_set_batch[year] = offshore
                
                imp_set_batch[year] = imp
                imp_nodal_set_batch[:, year] = imp_nodal
                h2imp_set_batch[year] = h2imp 
                h2imp_nodal_set_batch[:, year] = h2imp_nodal

                lly_set_batch[:, year] = lly
                plly_set_batch[:, year] = plly

                elinst_set_batch[year] = elinst
                elinst_nodal_set_batch[:, year] = elinst_nodal
                ellf_set_batch[year] = ellf
                ellf_nodal_set_batch[:, year] = ellf_nodal

                h2topinst_set_batch[year] = h2topinst 
                h2topinst_nodal_set_batch[:, year] = h2topinst_nodal
                h2toplf_set_batch[year] = h2toplf 
                h2toplf_nodal_set_batch[:, year] = h2toplf_nodal

                soc_h2_max_res_set_batch[year] = soc_h2_max_res  

                cost_set_batch[year] = cost
                runtimeopti_set_batch[year] = runtime_opti

            end # end year loop

## Thread lock
            Threads.lock(u) do

                #Time - everything locked
                start_blocktowrite = now()
                println("--- Thread ", Threads.threadid()," - Time : ",hour(start_blocktowrite), ":", minute(start_blocktowrite), ":", second(start_blocktowrite), " - Writing & all threads BLOCKED")

## Append batch variables to global variables

                global n_mc = n_mc + n_years_in_batch

                global lol_set = append!(lol_set, lol_set_batch)
                global lol_nodal_set = hcat(lol_nodal_set, lol_nodal_set_batch)
                global loe_set = append!(loe_set, loe_set_batch)
                global loe_nodal_set = hcat(loe_nodal_set, loe_nodal_set_batch)
                global h2loe_set = append!(h2loe_set, h2loe_set_batch)
                global h2loe_nodal_set = hcat(h2loe_nodal_set, h2loe_nodal_set_batch)

                global e_curt_set = append!(e_curt_set, e_curt_set_batch)
                global e_curt_nodal_set = hcat(e_curt_nodal_set, e_curt_nodal_set_batch)
                global p_curt_set = append!(p_curt_set, p_curt_set_batch)

                global phscycle_set = hcat(phscycle_set, phscycle_set_batch)
                global battcycle_set = hcat(battcycle_set, battcycle_set_batch)
                global h2cycle_set = hcat(h2cycle_set, h2cycle_set_batch)

                global shift_set = append!(shift_set, shift_set_batch)
                global shift_nodal_set = hcat(shift_nodal_set, shift_nodal_set_batch) 
                global shed_set = append!(shed_set, shed_set_batch)
                global shed_nodal_set = hcat(shed_nodal_set, shed_nodal_set_batch) 

                global conv_set = append!(conv_set, conv_set_batch)
                global genlf_set = hcat(genlf_set, genlf_set_batch)
                global offshore_set = append!(offshore_set, offshore_set_batch)
                
                global imp_set = append!(imp_set, imp_set_batch)
                global imp_nodal_set = hcat(imp_nodal_set, imp_nodal_set_batch)
                global h2imp_set = append!(h2imp_set, h2imp_set_batch)
                global h2imp_nodal_set = hcat(h2imp_nodal_set, h2imp_nodal_set_batch)

                global lly_set = hcat(lly_set, lly_set_batch)
                global plly_set = hcat(plly_set, plly_set_batch)

                global elinst_set = append!(elinst_set, elinst_set_batch)
                global elinst_nodal_set = hcat(elinst_nodal_set, elinst_nodal_set_batch)
                global ellf_set = append!(ellf_set, ellf_set_batch)
                global ellf_nodal_set = hcat(ellf_nodal_set, ellf_nodal_set_batch)

                global h2topinst_set = append!(h2topinst_set, h2topinst_set_batch)
                global h2topinst_nodal_set = hcat(h2topinst_nodal_set, h2topinst_nodal_set_batch)
                global h2toplf_set = append!(h2toplf_set, h2toplf_set_batch)
                global h2toplf_nodal_set = hcat(h2toplf_nodal_set, h2toplf_nodal_set_batch)

                global soc_h2_max_res_set = append!(soc_h2_max_res_set, soc_h2_max_res_set_batch)

                global cost_set = append!(cost_set, cost_set_batch)
                global runtimeopti_set = append!(runtimeopti_set, runtimeopti_set_batch)

## Compute expected values

                global lole = cumsum(lol_set) ./ (1:length(lol_set))
                global lole_nodal = cumsum(lol_nodal_set, dims=2) ./ transpose(repeat(1:size(lol_nodal_set)[2], 1, size(lol_nodal_set)[1]))
                global loee = cumsum(loe_set) ./ (1:length(loe_set))
                global loee_nodal = cumsum(loe_nodal_set, dims=2) ./ transpose(repeat(1:size(loe_nodal_set)[2], 1, size(loe_nodal_set)[1]))
                global h2loee = cumsum(h2loe_set) ./ (1:length(h2loe_set))
                global h2loee_nodal = cumsum(h2loe_nodal_set, dims=2) ./ transpose(repeat(1:size(h2loe_nodal_set)[2], 1, size(h2loe_nodal_set)[1]))

                global e_curte = cumsum(e_curt_set) ./ (1:length(e_curt_set))
                global e_curte_nodal = cumsum(e_curt_nodal_set, dims=2) ./ transpose(repeat(1:size(e_curt_nodal_set)[2], 1, size(e_curt_nodal_set)[1]))
                global p_curte = cumsum(p_curt_set) ./ (1:length(p_curt_set))

                global phscyclee = cumsum(phscycle_set, dims=2) ./ transpose(repeat(1:size(phscycle_set)[2], 1, size(phscycle_set)[1]))
                global battcyclee = cumsum(battcycle_set, dims=2) ./ transpose(repeat(1:size(battcycle_set)[2], 1, size(battcycle_set)[1]))
                global h2cyclee = cumsum(h2cycle_set, dims=2) ./ transpose(repeat(1:size(h2cycle_set)[2], 1, size(h2cycle_set)[1]))

                global shifte = cumsum(shift_set) ./ (1:length(shift_set))
                global shifte_nodal = cumsum(shift_nodal_set, dims=2) ./ transpose(repeat(1:size(shift_nodal_set)[2], 1, size(shift_nodal_set)[1])) 
                global shede = cumsum(shed_set) ./ (1:length(shed_set))
                global shede_nodal = cumsum(shed_nodal_set, dims=2) ./ transpose(repeat(1:size(shed_nodal_set)[2], 1, size(shed_nodal_set)[1]))

                global conve = cumsum(conv_set) ./ (1:length(conv_set))
                global genlfe = cumsum(genlf_set, dims=2) ./ transpose(repeat(1:size(genlf_set)[2], 1, size(genlf_set)[1]))
                global offshoree = cumsum(offshore_set) ./ (1:length(offshore_set))

                global impe = cumsum(imp_set) ./ (1:length(imp_set))
                global impe_nodal = cumsum(imp_nodal_set, dims=2) ./ transpose(repeat(1:size(imp_nodal_set)[2], 1, size(imp_nodal_set)[1]))
                global h2impe = cumsum(h2imp_set) ./ (1:length(h2imp_set))
                global h2impe_nodal = cumsum(h2imp_nodal_set, dims=2) ./ transpose(repeat(1:size(h2imp_nodal_set)[2], 1, size(h2imp_nodal_set)[1]))

                global llye = cumsum(lly_set, dims=2) ./ transpose(repeat(1:size(lly_set)[2], 1, size(lly_set)[1]))
                global pllye = cumsum(plly_set, dims=2) ./ transpose(repeat(1:size(plly_set)[2], 1, size(plly_set)[1]))

                global elinste = cumsum(elinst_set) ./ (1:length(elinst_set))
                global elinste_nodal = cumsum(elinst_nodal_set, dims=2) ./ transpose(repeat(1:size(elinst_nodal_set)[2], 1, size(elinst_nodal_set)[1])) 
                global ellfe = cumsum(ellf_set) ./ (1:length(ellf_set))
                global ellfe_nodal = cumsum(ellf_nodal_set, dims=2) ./ transpose(repeat(1:size(ellf_nodal_set)[2], 1, size(ellf_nodal_set)[1])) 

                global h2topinste = cumsum(h2topinst_set) ./ (1:length(h2topinst_set))
                global h2topinste_nodal = cumsum(h2topinst_nodal_set, dims=2) ./ transpose(repeat(1:size(h2topinst_nodal_set)[2], 1, size(h2topinst_nodal_set)[1]))
                global h2toplfe = cumsum(h2toplf_set) ./ (1:length(h2toplf_set))
                global h2toplfe_nodal = cumsum(h2toplf_nodal_set, dims=2) ./ transpose(repeat(1:size(h2toplf_nodal_set)[2], 1, size(h2toplf_nodal_set)[1])) 

                global soc_h2_max_rese = cumsum(soc_h2_max_res_set) ./ (1:length(soc_h2_max_res_set))

                global coste = cumsum(cost_set) ./ (1:length(cost_set))

                ## Update MC infos ##

                #global epsilon_loee_old = append!(epsilon_loee_old, last(epsilon_loee))

                if std(loee[1:end-1]) ==0 && mean(loee[1:end-1]) ==0
                    global epsilon_loee_old = append!(epsilon_loee_old, 0)
                else 
                    if convergence_type=="gamma"
                        global epsilon_loee_old = append!(epsilon_loee_old, std(loee[1:end-1])/(mean(loee[1:end-1]))) #gamma
                    elseif convergnce_type=="zeta"
                        global epsilon_loee_old = append!(epsilon_loee_old, std(loe_set[1:end-1])/(mean(loee[1:end-1])*sqrt(n_mc))) #zeta: attntion au n_mc... 
                    else
                        println("ERROR: CONVRGENCE TYPE NOT DEFINED")
                    end
                end

                if std(loee) ==0 && mean(loee) ==0
                    global epsilon_loee = append!(epsilon_loee, 0)
                else 
                    if convergence_type=="gamma"
                        global epsilon_loee = append!(epsilon_loee, std(loee)/(mean(loee))) #gamma
                    elseif convergnce_type=="zeta"
                        global epsilon_loee = append!(epsilon_loee, std(loe_set)/(mean(loee)*sqrt(n_mc))) #zeta
                    else
                        println("ERROR: CONVRGENCE TYPE NOT DEFINED")
                    end
                end

                if last(epsilon_loee) == 0 && last(epsilon_loee_old) == 0
                    global incr_epsilon_loee = append!(incr_epsilon_loee, 0)
                elseif last(epsilon_loee) != 0 && last(epsilon_loee_old) == 0
                    global incr_epsilon_loee = append!(incr_epsilon_loee,NaN)
                else
                    global incr_epsilon_loee = append!(incr_epsilon_loee, abs(last(epsilon_loee)-last(epsilon_loee_old))/last(epsilon_loee_old))
                end

                global incr_epsilon_loee_rm = rollingmedian(incr_epsilon_loee,30)

                println("epsilon_loee", last(epsilon_loee))
                println("epsilon_loee_old", last(epsilon_loee_old))
                println("incr_epsilon_loee", last(incr_epsilon_loee))
                println("incr_epsilon_loee_rm", last(incr_epsilon_loee_rm))

## Saving results in .txt files
                
                ## Write LOLE ## [h] 
                open("./Results/"*simu_name*"/output/lole.txt", "w") do io
                    writedlm(io, [round.(lol_set; digits=2) round.(lole; digits=2) ])   
                end

                ## Write nodal LOLE ## [h] 
                open("./Results/"*simu_name*"/output/lole_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(lol_nodal_set, lole_nodal)), digits=4))   
                end

                ## Write LOEE ## [GWh] 
                open("./Results/"*simu_name*"/output/loee.txt", "w") do io
                    writedlm(io, [round.(loe_set/1e3; digits=4) round.(loee/1e3; digits=4)])   
                end

                ## Write nodal LOEE ## [MWh] 
                open("./Results/"*simu_name*"/output/loee_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(loe_nodal_set, loee_nodal)), digits=4))   
                end

                ## Write H2 LOEE ## [h] 
                open("./Results/"*simu_name*"/output/h2loee.txt", "w") do io
                    writedlm(io, [round.(h2loe_set; digits=2) round.(h2loee; digits=2) ])   
                end

                ## Write nodal H2 LOEE ## [MWh_h2] 
                open("./Results/"*simu_name*"/output/h2loee_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(h2loe_nodal_set, h2loee_nodal)), digits=4))   
                end

                ## Write curtailment: E ## [GWh]
                open("./Results/"*simu_name*"/output/curtail.txt", "w") do io
                    writedlm(io, [round.(e_curt_set/1e3; digits=4) round.(e_curte/1e3; digits=4)])   
                end

                ## Write nodal curtailment: E ## [MWh]
                open("./Results/"*simu_name*"/output/curtail_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(e_curt_nodal_set, e_curte_nodal)), digits=4))   
                end

                ## Write peak curtailment: P ## [MW]
                open("./Results/"*simu_name*"/output/curtailpeak.txt", "w") do io
                    writedlm(io, [round.(p_curt_set; digits=4) round.(p_curte; digits=4)])   
                end

                ## Write phs infos ## [-]
                open("./Results/"*simu_name*"/output/phscycle.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(phscycle_set, phscyclee)), digits=4))   
                end

                ## Write batt infos ## [-]
                open("./Results/"*simu_name*"/output/battcycle.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(battcycle_set, battcyclee)), digits=4))   
                end

                ## Write h2 sto infos ## [-]
                open("./Results/"*simu_name*"/output/h2cycle.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(h2cycle_set, h2cyclee)), digits=4))   
                end

                ## Write shift infos ## [GWh]
                open("./Results/"*simu_name*"/output/shift.txt", "w") do io
                    writedlm(io, [round.(shift_set/1e3; digits=4) round.(shifte/1e3; digits=4)])   
                end

                ## Write nodal shift infos ## [MWh]
                open("./Results/"*simu_name*"/output/shift_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(shift_nodal_set, shifte_nodal)), digits=4))   
                end

                ## Write shed infos ## [GWh]
                open("./Results/"*simu_name*"/output/shed.txt", "w") do io
                    writedlm(io, [round.(shed_set/1e3; digits=4) round.(shede/1e3; digits=4)])   
                end

                ## Write nodal shed infos ## [MWh]
                open("./Results/"*simu_name*"/output/shed_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(shed_nodal_set, shede_nodal)), digits=4))   
                end

                ## Write conventional gen prod ## [GWh]
                open("./Results/"*simu_name*"/output/conv.txt", "w") do io
                    writedlm(io, [round.(conv_set/1e3; digits=4) round.(conve/1e3; digits=4)])   
                end

                ## Write conventional gen load factor ## [-]
                open("./Results/"*simu_name*"/output/genlf.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(genlf_set, genlfe)), digits=4))   #years = rows and cols=gens
                    #writedlm(io, [round.(genlf_set, digits=4) round.(genlfe, digits=4)]) #years=col and rows=gen  
                end

                ## Write offshore prod ## [GWh]
                open("./Results/"*simu_name*"/output/offshore.txt", "w") do io
                    writedlm(io, [round.(offshore_set/1e3; digits=4) round.(offshoree/1e3; digits=4)])   
                end

                ## Write imports ## [GWh]
                open("./Results/"*simu_name*"/output/imp.txt", "w") do io
                    writedlm(io, [round.(imp_set/1e3; digits=4) round.(impe/1e3; digits=4)])   
                end

                ## Write nodal imports  ## [MWh]
                open("./Results/"*simu_name*"/output/imp_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(imp_nodal_set, impe_nodal)), digits=4))   
                end

                ## Write h2 imports ## [GWh_h2]
                open("./Results/"*simu_name*"/output/h2imp.txt", "w") do io
                    writedlm(io, [round.(h2imp_set/1e3; digits=4) round.(h2impe/1e3; digits=4)])   
                end

                ## Write nodal h2 imports  ## [MWh_h2]
                open("./Results/"*simu_name*"/output/h2imp_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(h2imp_nodal_set, h2impe_nodal)), digits=4))   
                end

                ## Write nodal Yearly Line Loading  ## [-]
                open("./Results/"*simu_name*"/output/lly.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(lly_set, llye)), digits=4))   
                end

                ## Write nodal Yearly H2 PipeLine Loading  ## [-]
                open("./Results/"*simu_name*"/output/plly.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(plly_set, pllye)), digits=4))   
                end

                # Write electrolyzer installed infos ## [MW]
                open("./Results/"*simu_name*"/output/elinst.txt", "w") do io
                    writedlm(io, [round.(elinst_set; digits=4) round.(elinste; digits=4)])   
                end

                ## Write nodal electrolyzer installed infos ## [MW]
                open("./Results/"*simu_name*"/output/elinst_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(elinst_nodal_set, elinste_nodal)), digits=4))   
                end

                ## Write electrolyzer load factor ## [-]
                open("./Results/"*simu_name*"/output/ellf.txt", "w") do io
                    writedlm(io, [round.(ellf_set; digits=4) round.(ellfe; digits=4)])   
                end

                ## Write nodal electrolyzer load factor infos ## [-]
                open("./Results/"*simu_name*"/output/ellf_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(ellf_nodal_set, ellfe_nodal)), digits=4))   
                end

                ## Write hydrogen-to-power units installed infos ## [MW]
                open("./Results/"*simu_name*"/output/h2topinst.txt", "w") do io
                    writedlm(io, [round.(h2topinst_set; digits=4) round.(h2topinste; digits=4)])   
                end

                ## Write nodal hydrogen-to-power units installed infos ## [MW]
                open("./Results/"*simu_name*"/output/h2topinst_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(h2topinst_nodal_set, h2topinste_nodal)), digits=4))   
                end

                ## Write hydrogen-to-power units load factor ## [-]
                open("./Results/"*simu_name*"/output/h2toplf.txt", "w") do io
                    writedlm(io, [round.(h2toplf_set; digits=4) round.(h2toplfe; digits=4)])   
                end

                ## Write nodal hydrogen-to-power units load factor infos ## [-]
                open("./Results/"*simu_name*"/output/h2toplf_nodal.txt", "w") do io
                    writedlm(io, round.(transpose(vcat(h2toplf_nodal_set, h2toplfe_nodal)), digits=4))   
                end

                ## Write h2 storage installed infos ## [MWh_h2]
                open("./Results/"*simu_name*"/output/soc_h2_max_res.txt", "w") do io
                    writedlm(io, [round.(soc_h2_max_res_set; digits=4) round.(soc_h2_max_rese; digits=4)])   
                end

                ## Write cost ## [€]
                open("./Results/"*simu_name*"/output/cost.txt", "w") do io
                    writedlm(io, [round.(cost_set; digits=4) round.(coste; digits=4)])   
                end

                ## Write epsilon_loee - epsilon_loee_old - incr_epsilon_loee ## [-] 
                open("./Results/"*simu_name*"/output/epsilon_loee.txt", "w") do io
                    writedlm(io, [epsilon_loee epsilon_loee_old incr_epsilon_loee incr_epsilon_loee_rm])   
                end

                ## Write run times ## [€]
                open("./Results/"*simu_name*"/runtimeopti.txt", "w") do io
                    writedlm(io, runtimeopti_set )   
                end

                stop_simu = now()
                
                ## Write log ##
                open("./Results/"*simu_name*"/log.txt", "w") do io
                    write(io, "LOLE [h] : "*string(last(lole))) 
                    write(io, "\rLOEE [GWh] : "*string(last(loee)/1e3)) 
                    write(io, "\rNumber of MC years : "*string(n_mc))
                    write(io, "\rincr_epsilon_loee [-] : "*string(last(incr_epsilon_loee)))
                    write(io, "\rcost [€/year] : "*string(last(coste))) 
                    write(io, "\rincr_limit [-] : "*string(incr_limit)) 
                    write(io, "\rn_min [years] : "*string(n_min)) 
                    write(io, "\rn_max [years] : "*string(n_max))
                    write(io, "\rn_threads (number of tasks) : "*string(n_threads))
                    write(io, "\rn_years_in_batch : "*string(n_years_in_batch))
                    write(io, "\rnumber of --threads : "*string(Threads.nthreads())) 
                    write(io, "\relapsed_cal : "*string(canonicalize(Dates.CompoundPeriod(stop_simu-start_simu)))) 
                end

                stop_blocktowrite = now()
                runtime_blocktowrite = Dates.value(stop_blocktowrite - start_blocktowrite)/1000 #milliseconds to seconds
                println("--- Thread ", Threads.threadid()," - Time : ",hour(stop_blocktowrite), ":", minute(stop_blocktowrite), ":", second(stop_blocktowrite), " - Stop blocking & end writing : ", runtime_blocktowrite, " seconds")
            end # end thread lock

## Check convergence

            println("n_mc ", n_mc)
            if n_mc >= n_max || (last(incr_epsilon_loee_rm) <= incr_limit && n_mc >= n_min)
                return
            end

        end # end batch loop
    end # end thread loop
end # end function MC





