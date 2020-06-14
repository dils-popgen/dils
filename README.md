# DILS_web
# Table of contents
1. [Dependencies](#1---dependencies)     
2. [Installation](#2---installation)  
3. [Execution](#3---execution)  

# 1 - Dependencies
	[singularity](https://sylabs.io/docs/) (tested with 3.1.1) 
	
# 2 - Installation
	1. clone git repository  
	```
	git clone https://github.com/popgenomics/DILS_web  
	```

	2. move into the DILS repertory  
	```
	cd DILS_web  
	```
		
	3. build singularity image (could take ~20minutes):  
	```  
	sudo singularity build DILS.sif DILS.def  
	```

# 3 - Execution  	
	1. run the shiny app (may need root premissions):  
	```
	sudo singularity exec --bind DILS/:/mnt DILS.sif host=[ip adress of your server] port=[port number where shiny is reachable] nCPU=[maximum number of CPUs to use simultaneously]
	```    
	 eg with a big machine with 100 CPUs:  
	```  
	sudo singularity exec --bind DILS/:/mnt DILS.sif webinterface/app.R host=127.0.0.9 port=8912 nCPU=100
	```  
	
	Please keep in mind that the max number of CPUs is the maximum number of CPUs DILS will use at certain times, but that DILS will not use 100% of the indicated number of CPUs throughout its whole run. This maximum usage will be punctual.  
	  
	2. shiny app is now available in your web browser at http://[ip adress of your server]:[pourt number],  
	eg:  
	```
	http://127.0.0.9:8912/
	```

