# DILS
# Table of contents
1. [Dependencies](#1---dependencies)     
2. [Installation](#2---installation)  
3. [Execution](#3---execution)  
4. [Example](#4---example)  

# 1 - Dependencies
Only tested on different free Linux distributions, not payable OS.  
Singularity has to be installed on the machine.  
[singularity](https://sylabs.io/docs/) (tested with 3.1.1)   
	
# 2 - Installation
## clone git repository  
```
git clone https://github.com/popgenomics/DILS_web  
```

## move into the DILS repertory  
```
cd DILS_web  
```
	
## build singularity image (could take ~20minutes):  
```  
sudo singularity build DILS.sif DILS.def  
```

# 3 - Execution  	
## run the shiny app (may need root premissions):  
sudo singularity exec --bind DILS/:/mnt DILS.sif host=[ip adress of your server] port=[port number where shiny is reachable] nCPU=[maximum number of CPUs to use simultaneously]
  
eg with a big machine with 100 CPUs:  
```  
sudo singularity exec --bind DILS/:/mnt DILS.sif webinterface/app.R host=127.0.0.9 port=8912 nCPU=100
```  

Please keep in mind that the max number of CPUs is the maximum number of CPUs DILS will use at certain times, but that DILS will not use 100% of the indicated number of CPUs throughout its whole run. This maximum usage will be punctual.  
  
shiny app is now available in your web browser at http://[ip adress of your server]:[port number],  
eg:  
```
http://127.0.0.9:8912/
```
But chose the IP adress and port number you want 
  
# 4 - Example  
## Execute DILS in your web browser
In your terminal, from the DILS_web directory:  
```  
singularity exec --bind DILS/:/mnt DILS.sif webinterface/app.R host=127.0.0.9 port=8912 nCPU=100
```  
May require `sudo` to be run
  
Then in your web brower:  
```
http://127.0.0.9:8912/
```

## Running ABC from a fasta file  
All example files are in the **example** sub_directory.  

### 1. Unarchive the example file  
```
tar -Jxvf mytilus.tar.xz
```
Will generates the input file: 
```
mytilus.fas
```

###  2. Run a demographic ABC analysis  
Upload the **fasta** file by clicking on:  
1. ABC  
2. Upload data  
3. Browse (Input file upload)  
  
Then you can set up the analysis and execute it.  

## 3. Exploring the results
Upload the archive produced by DILS in a **tar.gz format** (doesn't need to be extracted) by clicking on:  
1. Results visualization  
2. Upload results  
3. Browse (Results to upload)  

Then you can explore the observed data and the results of the inferences.  

 
