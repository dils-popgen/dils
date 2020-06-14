# DILS_web

1. Dependencies
	- [singularity](https://sylabs.io/docs/) (tested with 3.1.1) 
	
2. Installation
	1. clone git repository  
	`git clone https://github.com/popgenomics/DILS_web`  

	2. move into the DILS repertory  
	`cd DILS_web`  
		
	3. build singularity image (could take ~20minutes):  
	`sudo singularity build DILS.sif DILS.def`  
	
	4. run the shiny app (may need root premissions):  
	`sudo singularity exec --bind DILS/:/mnt DILS.sif host=[ip adress of your server] port=[port number where shiny is reachable] nCPU=[maximum number of CPUs to use simultaneously]`  
	 eg:  
	`sudo singularity exec --bind DILS/:/mnt DILS.sif webinterface/app.R host=127.0.0.9 port=8912 nCPU=100`
	
	5. shiny app is now available in your web browser at http://[ip adress of your server]:[pourt number],  
	eg:  
	`http://127.0.0.9:8912/`

