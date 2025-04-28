
from neuron import h 
from neuron.units import ms, mV 
from matplotlib import pyplot as plt
import numpy as np
import math 
import pandas as pd 
from ecc_IB import PYRS 
from ini_IFS import INFS 
from RS_Cell import RS_CELL
h.load_file("stdrun.hoc")

def NetworkModel(sigma, alpha):
    
    N_tot = 100 # number of cells
    
    p_ecc = 0.8 # excitatory cells percentage 
    p_ini = 0.2 # inhibitory cells percentage 

    N_ecc = int(N_tot*p_ecc) # number of excitatory cells
    
    p_ecc_IB = 0.75 # percentage of Intrinsic Bursting 
    p_ecc_RS = 0.25 # percentage of Regular Spiking 
    N_ecc_IB = int(N_ecc * p_ecc_IB) 
    N_ecc_RS = int(N_ecc * p_ecc_RS)  
    
    N_ini = int(N_tot*p_ini) # Number of inhibitory neurons
    
    XEcc = np.random.rand(N_ecc) 
    YEcc = np.random.rand(N_ecc) 
    XIni = np.random.rand(N_ini) 
    YIni = np.random.rand(N_ini) 
    
    
    EE = np.zeros((N_ecc,N_ecc))
    EI = np.zeros((N_ecc,N_ini))
    II = np.zeros((N_ini,N_ini))
    IE = np.zeros((N_ini,N_ecc))
    
    EccCellsList = [] 

    for k in range(N_ecc_IB): 
        EccCellsList.append(PYRS())
    
    for k in range(N_ecc_IB,N_ecc,1): 
        EccCellsList.append(RS_CELL())    
    
    IniCellsList = [] 
    for k in range(N_ini): 
        IniCellsList.append(INFS()) 
    
    arr_ecc = np.arange(N_ecc) 
        
    arr_ini = np.arange(N_ini) 
        
    syn_E_E = [] 
    netcons_E_E = [] 
    
    syn_E_I = [] 
    netcons_E_I = [] 
    
    syn_I_E = [] 
    netcons_I_E = [] 
    
    syn_I_I = [] 
    netcons_I_I = [] 
    
    syn_mech_w_ecc = 0.025*0.25 
    syn_mech_w_RS = 0.025*0.25
    syn_mech_w_IB = 0.025*0.65
    syn_mech_w_ini = 0.02*0.25
    
    stim_dur = 61000 
    stim_start = 10
    T_stim = 100  

    p_EE = 0.1; p_EI = 0.3; p_II = 0.2; p_IE = 0.4; # connections probabilities 
    
    ampl=8/((np.sqrt(2)*np.sqrt(math.pi)*sigma/2)*math.erf(np.sqrt(2)/(2*sigma))*0.1*80)
    
    vreve=0
    vrevi=-70
    vt=-60
    epsp_ampl=0.7
    w_E_E=ampl*epsp_ampl/(vreve-vt) 
    w_E_I =w_E_E/3
    w_I_E=(-alpha*ampl*epsp_ampl/(vrevi-vt))
    w_I_I =w_I_E/3 
    
    
    w_RS = 0 

    for i in range(N_ecc): 
        
        xS = XEcc[i]; yS = YEcc[i] 
    
        new_ecc = np.delete(arr_ecc,np.where(arr_ecc == i)) 
    
        XNewEcc = XEcc[new_ecc] 
        YNewEcc = YEcc[new_ecc] 
        dist = np.sqrt(np.power((xS - XNewEcc),2) + np.power((yS - YNewEcc),2))
        PGaussEE = p_EE * np.exp(-( np.power((xS - XNewEcc),2) + np.power((yS - YNewEcc),2))/(2*np.power(sigma,2)) )
        PGaussEI = p_EI * np.exp(-( np.power((xS - XIni),2) + np.power((yS - YIni),2))/(2*np.power(sigma,2)))
        RndmFiltEE = np.random.rand(np.size(PGaussEE)) 
        RndmFiltEI = np.random.rand(np.size(PGaussEI)) 
    
        NewArrPercEE = new_ecc[RndmFiltEE<PGaussEE] 
        NewArrPercEI = arr_ini[RndmFiltEI<PGaussEI]
        EE[i,NewArrPercEE] = 1 
        EI[i,NewArrPercEI] = 1
    
        
        for j in NewArrPercEE:
            syn = h.Exp2SynM(EccCellsList[j].soma(0.5))
            syn.tau1 = 0.2
            syn.tau2 = 2.0
            syn.e = 0
            syn.dd = 0.8
            syn.taud = 10000
    
            nc_E_E = h.NetCon(EccCellsList[i].soma(0.5)._ref_v,syn,0,0,w_E_E,sec=EccCellsList[i].soma) 
            netcons_E_E.append(nc_E_E) 
            syn_E_E.append(syn)
    
    
        
        for j in NewArrPercEI:
            
            syn = h.Exp2SynM(IniCellsList[j].soma(0.5))
            syn.tau1 = 0.2 
            syn.tau2 = 2.0 
            syn.e = 0
            syn.dd = 0.8
            syn.taud = 10000
    
            nc = h.NetCon(EccCellsList[i].soma(0.5)._ref_v,syn,0,0,w_E_I,sec=EccCellsList[i].soma)
    
            netcons_E_I.append(nc)
            syn_E_I.append(syn) 
    
       
    spike_times_E_E = [h.Vector() for ncEE in netcons_E_E] 
    aus = zip(netcons_E_E, spike_times_E_E)
    del(xS,yS) 
    
    
    nc_st = [] 
    for i in range(N_ecc): 
        nc_aus = h.NetCon(EccCellsList[i].soma(0.5)._ref_v, None, sec=EccCellsList[i].soma)
        nc_aus.threshold = 0 
        nc_st.append(nc_aus)
    NewSpikeTimes = [h.Vector() for NCAus in nc_st]
    for NCAus, SpikeTimesVec in zip(nc_st,NewSpikeTimes): 
        NCAus.record(SpikeTimesVec)
    
    
    for i in range(N_ini): 
        
        xS = XIni[i]; yS = YIni[i] 
        new_ini = np.delete(arr_ini,np.where(arr_ini == i)) 
        XNewIni = XIni[new_ini] 
        YNewIni = YIni[new_ini] 
        PGaussII = p_II * np.exp(-( np.power((xS - XNewIni),2) + np.power((yS - YNewIni),2))/(2*np.power(sigma,2)) )
        PGaussIE = p_IE * np.exp(-( np.power((xS - XEcc),2) + np.power((yS - YEcc),2))/(2*np.power(sigma,2)) )
        RndmFiltII = np.random.rand(np.size(PGaussII))
        RndmFiltIE = np.random.rand(np.size(PGaussIE)) 
        NewArrPercII = new_ini[RndmFiltII<PGaussII]
        NewArrPercIE = arr_ecc[RndmFiltIE<PGaussIE]
        II[i,NewArrPercII] = 1 
        IE[i,NewArrPercIE] = 1 
    
        
        for j in NewArrPercII:
            
            syn = h.Exp2SynM(IniCellsList[j].soma(0.5))
            syn.tau1 =  0.1
            syn.tau2 = 3.0 
            syn.e = -70
            syn.dd = 0.8
            syn.taud = 10000
                
            nc_I_I = h.NetCon(IniCellsList[i].soma(0.5)._ref_v,syn,0,0,w_I_I,sec=IniCellsList[i].soma)
    
            netcons_I_I.append(nc_I_I)
            syn_I_I.append(syn)
    

        for j in NewArrPercIE:
            
            syn = h.Exp2SynM(EccCellsList[j].soma(0.5))
            syn.tau1 = 0.1 
            syn.tau2 = 3.0 
            syn.e =  -70
            syn.dd = 0.8
            syn.taud = 10000
            
            nc = h.NetCon(IniCellsList[i].soma(0.5)._ref_v,syn,0,0,w_I_E,sec=IniCellsList[i].soma)
    
            netcons_I_E.append(nc) 
            syn_I_E.append(syn)
    
         
    spike_times_I_I = [h.Vector() for ncII in netcons_I_I]
    
    nc_st_I = [] 
    for i in range(N_ini): 
        nc_aus_I = h.NetCon(IniCellsList[i].soma(0.5)._ref_v, None, sec=IniCellsList[i].soma) 
        nc_aus_I.threshold = 0 
        nc_st_I.append(nc_aus_I) 
        
    NewSpikeTimesIni = [h.Vector() for NCAusIni in nc_st_I]
    for NCAusIni, SpikeTimesVecIni in zip(nc_st_I,NewSpikeTimesIni): 
        NCAusIni.record(SpikeTimesVecIni)
    
    
    p_ecc_stim = 1 
    
    n_cell_ecc_stim = int(p_ecc_stim*N_ecc) 
                                       
    p_ini_stim = 1 
    n_cell_ini_stim = int(p_ini_stim*N_ini)
   
    arr_ecc_stim = np.random.choice(arr_ecc,n_cell_ecc_stim,replace=False)
    arr_ini_stim = np.random.choice(arr_ini,n_cell_ini_stim,replace=False)
    
    ''' 
    ::::::::::::::: Synaptic connections with NetStim :::::::::::::::::::::::::
    ''' 
    nu_stim = 1/ (T_stim*1e-3)
    
    arr_syn_stim_ecc = [] 
    arr_nc_netstim_ecc = [] 
    arr_stim_ecc = [] 
    for i in arr_ecc_stim: 
      
        if i < N_ecc_IB:
            syn_mech_w_ecc = syn_mech_w_IB 
        else: 
            syn_mech_w_ecc = syn_mech_w_RS 
            
        syn_stim_ecc = h.Exp2SynM(EccCellsList[i].soma(0.5))
        syn_stim_ecc.tau1 = 0.2
        syn_stim_ecc.tau2 = 2.0
        syn_stim_ecc.e = 0
        syn_stim_ecc.dd = 1.0
        
        stim_ecc = h.NetStim() 
        
        stim_ecc.number = (stim_dur-stim_start)/T_stim + 500/T_stim
        stim_ecc.start = stim_start
        stim_ecc.noise = 1
        stim_ecc.interval = T_stim 
       
        nc_netstim_ecc = h.NetCon(stim_ecc,syn_stim_ecc) 
        nc_netstim_ecc.delay = 0 
        nc_netstim_ecc.weight[0] = syn_mech_w_ecc
        
        arr_syn_stim_ecc.append(syn_stim_ecc)
        arr_nc_netstim_ecc.append(nc_netstim_ecc)
        
        arr_stim_ecc.append(stim_ecc)
        
    arr_syn_stim_ini = [] 
    arr_nc_netstim_ini = [] 
    arr_stim_ini = []; 
        
    for j in arr_ini_stim: 
        
        syn_stim_ini = h.Exp2SynM(IniCellsList[j].soma(0.5)) 
        syn_stim_ini.tau1 = 0.2
        syn_stim_ini.tau2 = 2.0
        syn_stim_ini.e = 0 
        syn_stim_ini.dd = 1.0
        
        stim_ini = h.NetStim() 
        
        stim_ini.number = (stim_dur-stim_start)/T_stim + 500/T_stim
        stim_ini.start = stim_start
        stim_ini.noise = 1 
        stim_ini.interval = T_stim 
    
        nc_netstim_ini = h.NetCon(stim_ini,syn_stim_ini) 
        nc_netstim_ini.delay = 0 
        nc_netstim_ini.weight[0] = syn_mech_w_ini 
        
        arr_syn_stim_ini.append(syn_stim_ini)
        arr_nc_netstim_ini.append(nc_netstim_ini)
        
        arr_stim_ini.append(stim_ini)
    v_vec_ecc = h.Vector() 
    v_vec_ini = h.Vector() 
    
    def RecV(cell_element): 
        rec_v = h.Vector() 
        rec_v.record(cell_element.soma(0.5)._ref_v)
        return rec_v
    

    h.step(0.025)
    h.finitialize(-70*mV)
    h.continuerun(stim_dur*ms)    
    
    spike_times_python = [] 
  
    
    for i, spike_times_vec in enumerate(spike_times_E_E):
       
        spike_times_python.append(list(spike_times_vec))
    NewSpikeTimesPy = [] 
    for i, SpikeTimesVec in enumerate(NewSpikeTimes): 
        NewSpikeTimesPy.append(list(SpikeTimesVec)) 
    
    spike_times_python_Ini = [] 
    for i, spike_times_vec in enumerate(spike_times_I_I): 
        spike_times_python_Ini.append(list(spike_times_vec)) 
    
    NewSpikeTimesIniPy = [] 
    for i, SpikeTimesVec in enumerate(NewSpikeTimesIni): 
        NewSpikeTimesIniPy.append(list(SpikeTimesVec)) 
        
        
    AllNewSpikeTimesPy = [] 
    for i in range(N_ecc): 
        AllNewSpikeTimesPy.append(NewSpikeTimesPy[i])
    for i in range(N_ini): 
        AllNewSpikeTimesPy.append(NewSpikeTimesIniPy[i])
    
    ID= 'alpha_' + str(alpha) + '_sigma_' + str(sigma)
    
    NomePlot = 'RasterPlot_' + str(N_tot) + '_' + ID + '.png'
    plt.figure(figsize=(14,6))
    plt.eventplot(AllNewSpikeTimesPy,linestyles='solid', 
               linewidths=3,color='blue')
    plt.xlabel('time (ms)',fontsize=30)
    plt.xticks(fontsize=30); plt.yticks(fontsize=30)
    plt.ylabel('Cell ID (#)',fontsize=30)
    plt.title('sigma = ' + str(sigma) +'; alpha = ' + str(alpha) + '; ampl = ' + str(ampl) + '\n' + 'syn_mec_w_ecc = ' + str(syn_mech_w_ecc) + '; syn_mech_w_ini = ' + str(syn_mech_w_ini) + '\n' + 'w_EE = ' + str(w_E_E) + '; w_EI = ' +
              str(w_E_I) + ';     w_IE =' + str(w_I_E) + '; w_II = ' + str(w_I_I))
    plt.xlim(2000,stim_dur)
    plt.savefig(NomePlot)
    plt.close()
    
    name1= 'AllNewSpikeTimes_' + str(N_tot) + '_' + ID + '.csv'
    DF2 = pd.DataFrame(AllNewSpikeTimesPy) 
    DF2.to_csv(name1)
    
