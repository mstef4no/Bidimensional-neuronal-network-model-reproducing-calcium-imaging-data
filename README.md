In silico neuronal cultures for Alzheimer's Disease study 

The codes allow to simulate the behavior of bidimensional neuronal cultures. 
A set of specific parameters can be varied to describe different scenarios according to the  specific research needs. 

**MainNetwork.py**

In the *MainNetwork.py* you can change the number of cells ```N_tot```;  
The relative weight of excitatory synapses to tune excitation/inhibition balance ```w_IE```, modifying the parameter ```alpha``` defined in the script *RunModel.py*;  
The spatial connectivity range defined in the script as ```P_EE```, ```P_EI```, ```P_II``` and ```P_IE```, modifying the parameter ```sigma``` in the script *RunModel.py*;  
The fraction of bursting cells, changing the percentage of Intrinsic Bursting cells ```p_ecc_IB``` and Regular Spiking cells "p_ecc_RS"; 
The contribution of NMDA synapses including or not the code sections related to the NMDA synaptic receptors contibution defined ion the excitatory cells section. 
The output of this script are the spike times of the network. 
