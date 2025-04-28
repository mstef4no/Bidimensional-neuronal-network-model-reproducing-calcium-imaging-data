
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
    
    
    #name2 = 'T_vec_' + str(N_tot)+ '_' + ID + '.csv'
    
    # cambio name1 e name2 
    #df2 = pd.read_csv(name2, dtype=np.float64)
    df1 = pd.read_csv(name1, dtype=np.float64) 
    
    df1_float = df1.astype(float) 
    #df2_float = df2.astype(float)
    #aus = df.to_numpy()
    #import numpy as geek 
    
    TempiDiSpike = np.array(df1) # , dtype= float)
    # arrotondare i dt di spike 
    # 
    
    #ArrayTemporale = np.array(df2)
    #ArrayTemporale = ArrayTemporale[0][1:]
    
    z = []
    dt_sim = 0.025 # delta t 
    T = 61000 # tempo finale in ms
    
    NewTimes = np.arange(0,T,0.025) # array temporale con step 0.025
    ArrayTemporale = NewTimes
    
    
    # Creo un treno di delta in corrispondenza dei tempi di spike, 
    # è necessario che i delta si trovino in corrispondenza degli stessi 
    # passi temporali fissati per la simulazione. Per questo si fa un arrotondamento 
    # dei tempi di spike 
    
    L = len(TempiDiSpike)
    for i in range(L): 
        aus = TempiDiSpike[i][1:] # in questo modo ho tolto la prima colonna che non serve 
        aus = aus[~(np.isnan(aus))]  # tolgo in NaN 
    
        # faccio un approssimazione dei tempi di spike in modo che rientrino nel passo temporale
        AA = aus / 0.025
        Ab=np.round(AA,0).astype(int)
    
        # creo un vettore di zeri 
        deltas_2 = np.zeros(int(T/0.025))
        
        # aggiungo gli 1 in corrispondenza dei tempi di spike. In questo modo il vettore ha la stessa grandezza 
        # di un generico vettore temporale della simulazione con passo 0.025 e durata 10000 ms 
        deltas_2[Ab]=1
        
        #plt.figure() 
        #plt.plot(deltas_2)
    
    
        z.append(deltas_2) # ha per ogni cellula un vettore della stessa lunghezza di quello temporale, 
        # ma con 0 e 1 in corrispondenza dei tempi di spike 
    
        del aus, AA, Ab, deltas_2 
    
    '''
    plt.figure() 
    plt.plot(t_calcio,z[5])
    '''
        
        
    # CALCOLO DELLA CONVOLUZIONE TRA IL KERNEL CALCIO E L'INSIEME DI DELTA 
    # IN CORRISPONDENZA DEI TEMPI DI SPIKE     
    tcaon=10; # 10 calcium rising time (milliseconds)
    tcaoff=500; # 700 calcium decay time (milliseconds)
    std_cnoise=0.02
    Aca=1;
    
    n_elem = T/dt_sim # numero di elementi 
    t_calcio = []
    t_calcio = np.arange(0,T,dt_sim) # 
    
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
            # da notare che facendo la convoluzione si ha l'unione dei vettori di
            # ascisse, quindi per avere l'effettiva convoluzione del kernel del calcio
            # seleziono solo la prima parte del vettore unione, fino alla lunghezza l
            
            #plt.figure()
            #plt.plot(ArrayTemporale[0],con[0:l])
            ResamplingCon = con[np.arange(i_cut,l,100)]
    
            noise_con = np.random.normal(0,std_cnoise,l-i_cut) #gaussian noise (mean,std,num_elements)
            noise_conR = np.random.normal(0,std_cnoise,len(ResamplingCon))
            
            ConvolArr.append(ResamplingCon+noise_conR) # insieme di vettori ricampionati
            NewConvolArr.append(con[i_cut:l]+noise_con) # insieme di vettori originali
        else:
            con = np.zeros((l-i_cut,),dtype=int)
            ResamplingCon = con[np.arange(0,l-i_cut,100)]
            
            noise_con = np.random.normal(0,std_cnoise,l-i_cut) #gaussian noise (mean,std,num_elements)
            noise_conR = np.random.normal(0,std_cnoise,len(ResamplingCon))
            
            ConvolArr.append(ResamplingCon+noise_conR) # insieme di vettori ricampionati
            NewConvolArr.append(con+noise_con) # insieme di vettori originali
    
        
        del con,ResamplingCon
        
    '''
    plt.figure() 
    plt.plot(ArrayTemporale[0], NewConvolArr[5])
    
    '''
        
        
        
        # manca quindi il calcolo della trasformata di fourier e dei coefficienti 
        # di correlazione come è stato fatto nei dati sperimentali 
    #a = range(1,5,1)
    #del a 
    #a = np.size(ArrayTemporale)
    #t_calcio_2 = t_calcio[np.arange(0,len(t_calcio),100)]
    
    ArrayTemporale2 = np.transpose(ArrayTemporale)
    
    # arr temporale della stessa lunghezza della convoluzione ricampionata 
    ArrayTemporale3 = ArrayTemporale2[np.arange(0,np.size(ArrayTemporale2),100)]
    ArrayTemporaleResamp = np.transpose(ArrayTemporale3)
    
    # DEVO CALCOLARE LA MEDIA NELLO SPAZIO 
    # quindi voglio un vettore che ha la stessa lunghezza di quello temporale 
    # ogni elemento è la media di quel valore su tutte le cellule 
    
    MeanConvol = np.nanmean(NewConvolArr, 0) # 0 indica la colonna, 1 la riga 
    
    MeanConvolResamp = np.nanmean(ConvolArr,0)
    
    
    '''
    plt.figure()
    plt.plot(ArrayTemporale3[i_cut_res:], np.transpose(ConvolArr))
    plt.xlim(t_cut,T)
    plt.show()
    
    plt.figure()
    plt.plot(ArrayTemporale[i_cut:], MeanConvol)
    plt.xlim(2000,T)
    plt.show()
    '''
    
    
    nameFig= 'MeanConvolRes_' + str(N_tot) + '_' + ID + '.png'
    plt.figure() 
    plt.plot(ArrayTemporale3[i_cut_res:], MeanConvolResamp)
    plt.title('Mean Convolution')
    plt.xlim(1000,T)
    plt.savefig(nameFig)
    plt.close()
    #plt.show()
    
    
    mean = np.mean(MeanConvol,0)
    NewMeanConvol = MeanConvol - mean 
    
    meanResamp = np.mean(MeanConvolResamp,0) 
    NewMeanConvolResamp = MeanConvolResamp - meanResamp 
    #DF = pd.DataFrame(NewMeanConvol)
    #DF.to_csv("NewMeanConvol.csv") 
    
    
    #np.savetxt('newconv.txt', NewConvolArr)
    #np.savetxt('newconv2.txt', ConvolArr)
    #nomeconv = 'Caconv_' + ID + '.txt'
    #np.savetxt(nomeconv,(ArrayTemporale[i_cut:],MeanConvol))
    nomeconv2 = 'CaconvRes_' + ID + '.txt'
    np.savetxt(nomeconv2,(ArrayTemporale3[i_cut_res:],MeanConvolResamp))

    
    # ----- :::::  ------ ------- -------- ------- ------- 
    Res = np.corrcoef(NewConvolArr) # matrice di coefficienti, completa 
    
    # ConvolArr è quella con i ricampionamento
    ResResamp = np.corrcoef(ConvolArr)
    
    UpTrResResamp = np.triu(ResResamp,k=1)
    UpTrResResamp[np.tril_indices(UpTrResResamp.shape[0], 0)] = np.nan
    UpTrRes = np.triu(Res,k=1)
    UpTrRes[np.tril_indices(UpTrRes.shape[0], 0)] = np.nan # mette dei nan 
    # nella triangolare inferiore, lo 0 indica da quale inizia a prendere la diagonale 
    
    BinSpacing = 0.025 
    BinEdges = np.arange(-1,1+BinSpacing,BinSpacing)
    # hist = np.histogram(UpTrRes)
    
    
    ResArr = [] 
    Res_old = Res 
    Res = []; 
    
    #Res = UpTrRes; # Devo considerare la triangolare senza gli 1
    Res = UpTrResResamp; # Devo considerare la triangolare senza gli 1
    
    #for i in range(Res.shape[0]): 
    #a_aus = Res[i] 
    ResArr = [Res[i] for i in range(Res.shape[0])]
    ResArrConc = np.concatenate((ResArr), axis=None)
    # a.append(aa)
    # 
    NomePlot = 'Istogramma_' + ID + '.png' 
    
    #ResArrConc_ = ResArrConc/max(ResArrConc)
    
    plt.figure(figsize=(8,5))
    R = plt.hist(ResArrConc, bins=BinEdges,color='blue') #, density=True, align='right') #, cumulative=False) 
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(NomePlot)
    plt.close()
    #plt.show()
   
    nomehist = 'CorrHist_' + ID + '.txt'
    np.savetxt(nomehist,ResArrConc)
    nomehist2 = 'CorrHistEdge_' + ID + '.txt'
    np.savetxt(nomehist2,BinEdges)

   #R = plt.hist(ResArrConc, bins=BinEdges,density=True)
    #RNorm = R[0]/max(R[0])
    #height = R[1]
    #RNorml = np.append(RNorm, 0)
    #plt.figure()
    #plt.bar(height, RNorml)
    
    #plt.ylim(0, 35)
    #plt.xlim(-0.5, -0.2) 
    #plt.savefig(NomePlot)
    
