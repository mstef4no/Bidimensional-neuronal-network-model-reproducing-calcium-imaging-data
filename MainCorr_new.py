
from scipy.signal import unit_impulse, convolve 
import numpy as np 
import pandas as pd 
import csv 
import matplotlib.pyplot as plt 

def NetworkCorr(sigma, alpha):
    
    IDload='alpha_' + str(alpha) + '_sigma_' + str(sigma)
    
    N_tot = 100
    name1= 'AllNewSpikeTimes_' + str(N_tot) + '_' + IDload + '.csv'
    
    ID='alpha_' + str("{:05.2f}".format(alpha)) + '_sigma_' + str("{:05.2f}".format(sigma))

    df1 = pd.read_csv(name1, dtype=np.float64) 
    
    df1_float = df1.astype(float) 
    
    TempiDiSpike = np.array(df1) 
    
    z = []
    dt_sim = 0.025 # delta t 
    T = 61000 # tempo finale in ms
    
    NewTimes = np.arange(0,T,0.025) 
    ArrayTemporale = NewTimes
    
    L = len(TempiDiSpike)
    for i in range(L): 
        aus = TempiDiSpike[i][1:] 
        aus = aus[~(np.isnan(aus))]  

        AA = aus / 0.025
        Ab=np.round(AA,0).astype(int)
    

        deltas_2 = np.zeros(int(T/0.025))
        
        deltas_2[Ab]=1
        
        z.append(deltas_2) 
        
        del aus, AA, Ab, deltas_2 
    
    tcaon=10; 
    tcaoff=500; 
    std_cnoise=0.02
    Aca=1;
    
    n_elem = T/dt_sim 
    t_calcio = []
    t_calcio = np.arange(0,T,dt_sim) 
    
    l = np.size(ArrayTemporale)
    ConvolArr = []; NewConvolArr = []; ResamplingCon = [];
    
    t_cut=1000
    i_cut=int(t_cut/dt_sim)
    i_cut_res=int(t_cut/(dt_sim*100))
    
    for i in range(L): 
    
        if np.sum(z[i])>0:
            y=(Aca*np.exp(-(t_calcio)/tcaoff))*(1-np.exp(-(t_calcio)/tcaon));
            con = convolve(z[i], y, mode='full') #, method='direct')
            con = con/np.max(con)
            
            ResamplingCon = con[np.arange(i_cut,l,100)]
    
            noise_con = np.random.normal(0,std_cnoise,l-i_cut) 
            noise_conR = np.random.normal(0,std_cnoise,len(ResamplingCon))
            
            ConvolArr.append(ResamplingCon+noise_conR) 
            NewConvolArr.append(con[i_cut:l]+noise_con) 
        else:
            con = np.zeros((l-i_cut,),dtype=int)
            ResamplingCon = con[np.arange(0,l-i_cut,100)]
            
            noise_con = np.random.normal(0,std_cnoise,l-i_cut) 
            noise_conR = np.random.normal(0,std_cnoise,len(ResamplingCon))
            
            ConvolArr.append(ResamplingCon+noise_conR) 
            NewConvolArr.append(con+noise_con) 
    
        
        del con,ResamplingCon
   
    ArrayTemporale2 = np.transpose(ArrayTemporale)
    
   
    ArrayTemporale3 = ArrayTemporale2[np.arange(0,np.size(ArrayTemporale2),100)]
    ArrayTemporaleResamp = np.transpose(ArrayTemporale3)
    
    MeanConvol = np.nanmean(NewConvolArr, 0) 
    
    MeanConvolResamp = np.nanmean(ConvolArr,0)
    
    nameFig= 'MeanConvolRes_' + str(N_tot) + '_' + ID + '.png'
    plt.figure() 
    plt.plot(ArrayTemporale3[i_cut_res:], MeanConvolResamp)
    plt.title('Mean Convolution')
    plt.xlim(1000,T)
    plt.savefig(nameFig)
    plt.close()
    
    mean = np.mean(MeanConvol,0)
    NewMeanConvol = MeanConvol - mean 
    
    meanResamp = np.mean(MeanConvolResamp,0) 
    NewMeanConvolResamp = MeanConvolResamp - meanResamp 
    
    nomeconv2 = 'CaconvRes_' + ID + '.txt'
    np.savetxt(nomeconv2,(ArrayTemporale3[i_cut_res:],MeanConvolResamp))

    
    
    Res = np.corrcoef(NewConvolArr) 
    

    ResResamp = np.corrcoef(ConvolArr)
    
    UpTrResResamp = np.triu(ResResamp,k=1)
    UpTrResResamp[np.tril_indices(UpTrResResamp.shape[0], 0)] = np.nan
    UpTrRes = np.triu(Res,k=1)
    UpTrRes[np.tril_indices(UpTrRes.shape[0], 0)] = np.nan 
    
    BinSpacing = 0.025 
    BinEdges = np.arange(-1,1+BinSpacing,BinSpacing)
    
    ResArr = [] 
    Res_old = Res 
    Res = []; 
    
    
    Res = UpTrResResamp; 
    
    
    ResArr = [Res[i] for i in range(Res.shape[0])]
    ResArrConc = np.concatenate((ResArr), axis=None)
    
    NomePlot = 'Istogramma_' + ID + '.png' 
    
    
    plt.figure(figsize=(8,5))
    R = plt.hist(ResArrConc, bins=BinEdges,color='blue') 
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(NomePlot)
    plt.close()
    
   
    nomehist = 'CorrHist_' + ID + '.txt'
    np.savetxt(nomehist,ResArrConc)
    nomehist2 = 'CorrHistEdge_' + ID + '.txt'
    np.savetxt(nomehist2,BinEdges)

    
