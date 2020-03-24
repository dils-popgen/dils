# DILS_web

1. Dependencies
	- [singularity](https://sylabs.io/docs/) (tested with 3.1.1) 
	- R (tested with 3.6.2)
		- R dependencies:
		- conda install  -c conda-forge r-shinycssloaders
		- conda install  -c conda-forge r-shinythemes
		- conda install  -c conda-forge r-shinydashboard
		- conda install  -c conda-forge r-shinydashboardplus
        - conda install  -c conda-forge r-shinyjs
        - conda install  -c conda-forge r-dt
        - conda install  -c conda-forge r-shinywidgets
        - conda install  -c conda-forge r-devtools
        - Rscript -e 'library(devtools); install_github("nik01010/dashboardthemes")'
        - conda install  -c conda-forge r-shinyhelper
        - conda install  -c conda-forge r-tidyr
        - conda install  -c conda-forge r-rcolorbrewer
        - conda install  -c r r-yaml
        - conda install  -c r r-data.table
        - conda install  -c conda-forge r-factominer
        - conda install  -c conda-forge r-ggplot2=3.2.1
        - conda install  -c conda-forge r-ggpubr
        - conda install  -c r r-nnet
        - conda install  -c conda-forge r-plotly
        - conda install  r-tidyverse
        - conda install  -c conda-forge r-viridis
        - conda install  -c bioconda r-matrixstats
        - conda install  -c conda-forge r-ranger
        - conda install  -c conda-forge r-rcpparmadillo
        - Rscript -e "install.packages('abcrf', repos='https://cloud.r-project.org/')"
	
2. Installation
	1. clone git repository
	2. build singularity image: `singularity build popgenomics.sif popgenomics.def`
	3. run the shiny app: `./popgenomics.sif host=[ip adress of your server] port=[port number where shiny is reachable]` eg: `./popgenomics.sif host=127.0.0.1 port=500`
	4. shiny app is now available in your web browser at http://[ip adress of your server]:[pourt number]
