import numpy as np
 
sigma_arr = np.array([0.1, 0.2,  0.5, 1.0, 10.0]) 
alpha_arr = np.array([0.01, 0.1, 0.15, 0.2, 1.0]) 



for i in range(len(sigma_arr)): 
    for j in range(len(alpha_arr)): 
    
        sigma = sigma_arr[i]; 
        alpha = alpha_arr[j]
        print(sigma, alpha)
        exec(open("./MainCorr.py").read()) 
        NetworkCorr(sigma, alpha)
