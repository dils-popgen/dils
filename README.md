# DILS_web

1. Dependencies
	- [singularity](https://sylabs.io/docs/) (tested on 3.1.1) 
2. Installation
	1. clone git repository
	2. build singularity image: `singularity build popgenomics.sif popgenomics.def`
	3. run the shiny app: `./popgenomics.sif host=[ip adress of your server] port=[port number where shiny is reachable]` 

eg: `./popgenomics.sif host=127.0.0.1 port=500`
	4. shiny app is now available in your web browser at http://[ip adress of your server]:[pourt number]
