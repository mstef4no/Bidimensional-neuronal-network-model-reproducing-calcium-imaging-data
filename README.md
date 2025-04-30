

The codes allow to simulate the behavior of bidimensional neuronal cultures. 
A set of specific parameters can be varied to describe different scenarios according to the  specific research needs. 

**MainNetwork.py**

In the *MainNetwork.py* you can change the number of cells ```N_tot```;  
The relative weight of excitatory synapses to tune excitation/inhibition balance ```w_IE```, modifying the parameter ```alpha``` defined in the script *RunModel.py*;  
The spatial connectivity range defined in the script as ```P_EE```, ```P_EI```, ```P_II``` and ```P_IE```, modifying the parameter ```sigma``` in the script *RunModel.py*;  
The fraction of bursting cells, changing the percentage of Intrinsic Bursting cells ```p_ecc_IB``` and Regular Spiking cells ```p_ecc_RS```; 
The contribution of NMDA synapses, including or not the code sections related to the NMDA synaptic receptors contribution defined in the excitatory cells section. 
The **output** of this script is the spike times of the network. 

**RunModel.py**
The "RunModel.py" is used to manage the "MainNetwork.py". 
In that script, it is possible to vary the ```alpha``` and ```sigma``` parameters to be used in the *MainNetwork.py* script. 

**RunCorr.py**  
The *RunCorr.py* script is used to manage the *MainCorr.py*.  
That script is used to extract and process the data provided by the *NetworkModel.py*.  
The **input** of that script are the spike times obtained from the *NetworkModel.py*.  
The **output** of that script are the calcium signals reconstruction, the Probability Spectral Analysis and the correlation coefficients.

