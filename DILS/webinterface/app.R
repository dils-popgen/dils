#!/usr/bin/Rscript

print(getwd())
#################################################################################################################################
#################################################################################################################################
#####                                                                                                                       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.                                          #####
#####                                                                                                                       #####   
#####    DILS is free software: you can redistribute it and/or modify                                                       #####
#####    it under the terms of the GNU General Public License as published by                                               #####
#####    the Free Software Foundation, either version 3 of the License, or                                                  #####
#####    (at your option) any later version.                                                                                #####
#####                                                                                                                       #####    
#####    DILS is distributed in the hope that it will be useful,                                                            #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                      #####
#####    GNU General Public License for more details.                                                                       #####
#####                                                                                                                       #####    
#####    You should have received a copy of the GNU General Public License                                                  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.                                                     #####
#####                                                                                                                       #####    
#####    Please send bugreports with examples or suggestions to                                                             #####
#####    camille.roux@univ-lille.fr                                                                                         #####
#####                                                                                                                       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####                                                                                                                       #####
#################################################################################################################################
#################################################################################################################################

for(tmp in commandArgs()){
	tmp = strsplit(tmp, '=')
	if(tmp[[1]][1] == 'port'){ port = as.numeric(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'host'){ host = as.character(tmp[[1]][2]) }
}

#options('shiny.port'=as.numeric(commandArgs()[2]), 'shiny.host'=commandArgs()[3])

library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(shinyWidgets)
library(dashboardthemes) # library(devtools); install_github("nik01010/dashboardthemes")
library(shinyhelper)
library(plotly)
library(viridis)
library(tidyr)
library(RColorBrewer)
library(yaml)
library(ggpubr)
library(FactoMineR)

pvalue = function(distribution, obs){
	obs = as.numeric(obs)
	distribution = as.numeric(distribution)
	median_x = median(distribution)
	if(obs==0 && median_x==0){
		return(NA)
	}else{
		if(obs==median_x){
			pval=0.5
		}else if(as.numeric(obs)>median_x){
			pval = length(which(distribution>as.numeric(obs)))/length(distribution)
		}else{
			pval = length(which(distribution<as.numeric(obs)))/length(distribution)
		}
	}
	return(pval)
}

plot3var = function (x, y, z, xlab = "", ylab = "", zlab = "", main = "", cex.lab = 1, couleurs = c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"), zlim = NULL, watermark = F, nlevels = 10, FUN="median"){
	median_z = c()
	mat = matrix(NA, length(table(y)), length(table(x)))
	colnames(mat) = names(table(x))
	rownames(mat) = names(table(y))
	ligne = 0
	colonne = 0
	for (x_i in as.numeric(names(table(x)))) {
		colonne = colonne + 1
		ligne = 0
		for (y_i in as.numeric(names(table(y)))) {
			ligne = ligne + 1
			if(FUN=='median'){
				mat[ligne, colonne] = median(z[which(x == x_i & y == y_i)])
				}
			if(FUN=='mean'){
				mat[ligne, colonne] = mean(z[which(x == x_i & y == y_i)])
			}
			if(FUN=='max'){
				mat[ligne, colonne] = max(z[which(x == x_i & y == y_i)])
			}
			if(FUN=='min'){
				mat[ligne, colonne] = min(z[which(x == x_i & y == y_i)])
			}
			if(FUN=='sd'){
				mat[ligne, colonne] = sd(z[which(x == x_i & y == y_i)])
			}
			median_z = c(median_z, mat[ligne, colonne])
		}
	}
	min_arr = which(mat == min(mat), arr.ind = T)
	max_arr = which(mat == max(mat), arr.ind = T)
	min_x = min_arr[, 2]
	max_x = max_arr[, 2]
	min_y = min_arr[, 1]
	max_y = max_arr[, 1]
	min_z = min(mat)
	max_z = max(mat)
	gradient = colorRampPalette(couleurs)
	dev.new(width = 8, height = 7)
	layout(matrix(c(1, 2), byrow = T, ncol = 2), width = c(4/5, 1/5))
	par(mar = c(4.5, 4, 4, 1), las = 1)
	if (is.null(zlim)) {
		zlim = range(mat)
	}

	if(length(min_x) == 1 && length(max_x) == 1){
		image(t(mat), xlab = "", ylab = "", col = gradient(nlevels), 
		cex.axis = cex.lab, axes = F, zlim = zlim)
		mtext(side = 3, text = main, line = 0.75, cex = cex.lab)
	}else{
		image(t(mat), xlab = "", ylab = "", col = gradient(nlevels), 
		cex.axis = cex.lab, axes = F, zlim = zlim)
		mtext(side = 3, text = main, line = 0.75, cex = cex.lab)
	}

	if (is.null(colnames(mat))) {
		mtext(side = 1, text = xlab, line = 2.5, cex = cex.lab)
		par(las = 3)
		mtext(side = 2, text = ylab, line = 2.75, cex = cex.lab)
	}
	else {
		migRates = rownames(mat)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at = 0:(length(table(x)) - 1)/(length(table(x)) - 1), labels = names(table(x)))
		mtext(xlab, 1, line = 2.5, cex = cex.lab)
		extRates = colnames(mat)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at = 0:(length(table(y)) - 1)/(length(table(y)) - 1), labels = names(table(y)))
		par(las = 0)
		mtext(ylab, 2, line = 2.75, cex = cex.lab)
	}
	if (watermark) {
		watermark()
	}
	
	par(las = 1)
	image.scale(mat, horiz = F, col = gradient(nlevels), xlab = "", ylab = "", cex.lab = cex.lab, cex.axis = cex.lab, zlim = zlim)
	par(las = 3)
	mtext(side = 2, text = zlab, line = 2.5, cex = cex.lab)
}


convertMenuItem <- function(mi,tabName) {
	mi$children[[1]]$attribs['data-toggle']="tab"
	mi$children[[1]]$attribs['data-value'] = tabName
	if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
		mi$attribs$class=NULL
	}
	mi
}

# welcome
welcome_page <- fluidPage(
	fluidRow(
		boxPlus(title = h2("Overview"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
			h3(strong("DILS"), "= ", strong("D", .noWS='outside'),"emographic ", strong("I", .noWS='outside'), "nferences with ", strong("L", .noWS='outside'), "inked ", strong("S", .noWS='outside'), "election"),
			h3(strong("DILS"), "is a DNA sequence analysis workflow to study the demographic history of sampled populations or species by using Approximate Bayesian Computations."),
			h3("From a single uploaded input file containing sequenced genes or DNA fragments,", strong("DILS"), "will:"),
			HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>1.</b> simulate different models/scenarios</h3>'),
			HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>2.</b> select the best model using an ABC approach based on <a href="https://cran.r-project.org/web/packages/abcrf/index.html" target="_blank"><font color="#c7f464"><b>random forests</b></font></a></h3>'),
			HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>3.</b> estimate the parameters of the best model using a <a href="https://cran.r-project.org/web/packages/abc/index.html" target="_blank"><font color="#c7f464"><b>neural network</b></font></a> approach</h3>'),
			HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>4.</b>measure the robustness of the analyses <b>DILS</b> is transparent on the ability of its inferences to reproduce the observed data or not</h3>'),
			hr(),
			h3("The primary goal of", strong("DILS"), "is to distinguish between isolation versus migration models of divergence between sister gene pools."),
			h3("Its ultimate goal is to produce for each studied gene the probability of being associated with a species barrier."),
			hr(),
			h3("Users with sequences data for a single group of individuals can also investigate alternative models of demographic change by using", strong(" DILS", .noWS='outside'), ".")
		)
	),
	
	fluidRow(
		#box(title = h2("Compared demographic models"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
		boxPlus(title = h2("How to use DILS?"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
			mainPanel(width=NULL,
				HTML('<h3>DILS is organized into two distinct features that can be seen on the side menu bar: <font color="#c7f464"><b>ABC</b></font> and <font color="#c7f464"><b>Results visualization</b></font>.</h3>'),
				htmlOutput('overview_DILS_picture'),
				HTML('<h3>The <b>ABC</b> feature comprises five tabs:</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Four to configure the ABC analysis (<i class="fas fa-cloud-upload-alt"></i>Upload Data, <i class="fas fa-bath"></i>Data Filtering, <i class="fas fa-users-cog"></i>Populations/Species and <i class="fas fa-dice"></i>Prior Distributions)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The last to execute from one to five ABC analyses (<i class="fas fa-microchip"></i>Run ABC)</h3>'),
				
				HTML('<h3>DILS requires a single input file whose format details are given in <i class="fas fa-industry"></i> / <i class="fas fa-cloud-upload-alt"></i></h3>'),
				HTML('<h3>In order to run the ABC analysis, it is necessary to validate your choices of configurations by clicking on the <font color="#c7f464"><b>Please check/validate your choices</b></font> button at the bottom of the four configuration pages</h3>'),
				HTML('<br><h3>The <b>Results visualization</b> feature comprises three tabs:</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To upload and preview the archive produced by DILS after the ABC analysis (<i class="fas fa-cloud-upload-alt"></i>Upload results)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To explore the <b>observed summary statistics</b> from the genomic data and the results of the <b>demographic inferences</b> obtained by the ABC (<i class="fas fa-lock"></i>User&#39s dataset)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To compare the results with previous analysis (<i class="fas fa-lock-open"></i>Collaborative science)</h3>')
				
			)
		)
	),
	
	fluidRow(
		#box(title = h2("Compared demographic models"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
		boxPlus(title = h2("Compared demographic models"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
			mainPanel(width=NULL, 
				h3(strong("1 population/species")),
				htmlOutput("models_picture_1pop"),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Constant</b> = single panmictic population of effective size <b><i>Ne</i></b> constant over time.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Expansion</b> = the size of the current population has suddenly become larger than in the past <b><i>T<sub>dem</sub></i></b> generations ago.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Contraction</b> = the current population experienced a decline in its size <b><i>T<sub>dem</sub></i></b> generations ago.</h3>'),
				hr(),
				h3(strong("2 populations/species")),
				htmlOutput("models_picture_2pop"),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>SI</b> = strict isolation: subdivision of an ancestral diploid panmictic population (of size <i>N<sub>anc</sub></i>) in two diploid populations (of constant sizes <i>N<sub>pop1</sub></i> and <i>N<sub>pop2</sub></i>) at time <i>T<sub>split</sub></i>.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>AM</b> = ancestral migration: the two newly formed populations continue to exchange alleles until time T<sub>AM</sub>.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>IM</b> = isolation with migration: the two daughter populations continuously exchange alleles until present time.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>SC</b> = secondary contact: the daughter populations first evolve in isolation (forward in time), then experience a secondary contact and start exchanging alleles at time T<sub>SC</sub>. Red phylogenies represent possible gene trees under each alternative model.</h3>'),
				hr(),
				h3(strong("4 populations/species")),
				htmlOutput("models_picture_4pop"),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A single generalist model, declined in 64 sub-models according to if there is one:</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;migration (bidirectional) or no migration between A and B, and/or between C and D.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bidirectional, unidirectional or no migration between A and C, and/or between B and D.</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In case of migration between A and B (and between C and D), the gene flow takes place since their separation T<sub>split_AB</sub> (and T<sub>split_CD</sub>).</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In case of migration between A and C (and between B and D), the gene flow occurs during a secondary contact T<sub>SC_AC</sub> (and T<sub>SC_BD</sub>) lower than min(c(T<sub>split_AB</sub>, T<sub>split_CD</sub>))</h3>')
			)
	),
	
	boxPlus(title = h2("Compared genomic models"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
		tags$style(type = "text/css", HTML(".irs-single { color: #1e2b37; font-size: 18px; background: transparent }")), # change the font of the sliderInput
		tags$head(tags$style(type='text/css', ".irs-grid-text { color: #556270; font-size: 12pt; }")),
		
		fluidRow(
			column(width=12,
				h3("All demographic models exist under AT LEAST two alternative genomic models concerning the effective size:"),
				HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) genomic homogeneity</b>, where the effective size <b><i>Ne</i></b> is genomically homogeneous (purple bar), <i>i.e.</i>, all locus are simulated by sharing the same <b>Ne</b> value. In this model <b>DILS</b> will try to estimate the value of <i>Ne</i> best explaining the observed data. <i>Ne</i> being independently estimated in all populations (current, past).'),
				HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) genomic heterogeneity</b> where <i>Ne</i> is genomically heterogeneous (yellow distribution), <i>i.e.</i>, all locus are simulated with a value of <i>Ne</i> drawn in a Beta(&alpha;, &beta;) distribution. In this model <b>DILS</b> will try to estimate the value of <i>Ne</i> as well as the two shape parameters <b>shape1</b> and <b>shape2</b> (&alpha; and &beta;) that best explain the observations. Here, <b>DILS</b> assumes that all populations (current and past) share the same Beta(&alpha;, &beta;) distribution but are independently rescaled by different <i>Ne</i> values.')
			),
			
			# Sidebar panel for inputs
			column(width=4,
				# Input: Slider for the number of bins
				chooseSliderSkin("HTML5", color = "#c7f464"),
			#	setSliderColor(rep("#556270",10), 1:10),
				sliderInput(inputId = "Ne", label = HTML('<h3>Effective population size:</h3>'), min = 0, max = 1000000, value = 10000, step=1000),
				sliderInput(inputId = "alpha", label = HTML('<h3>Shape parameter &alpha;:</h3>'), min = 1, max = 20, value = 10, step=0.1),
				sliderInput(inputId = "beta", label = HTML('<h3>Shape parameter &beta;:</h3>'), min = 1, max = 20, value = 3, step=0.1)	
			),
			
			# Main panel for displaying outputs
			column(width=8,
				# Output: density
				plotOutput(outputId = "genomic_hetero")
			)
		),
	
		fluidRow(
			column(width=12,
				hr(),
				h3("In addition, all demographic models with migration have two alternative models of introgression:"),
				HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) genomic homogeneity, </b>where all loci share the same introgression rate for a given direction. Here, <b>DILS</b> will simply try to <b>independently</b> estimate the introgression rate of each direction.</h3>'),
				HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) genomic heterogeneity, </b>where introgression rates vary throughout the genome. This variation can follow a <b>beta</b> or <b>bimodal</b> distribution (see below).</h3>'),
				htmlOutput("homo_hetero")
	
			)
		),
		
		fluidRow(
			column(
				hr(),
				width = 12,
				HTML('<h3>Genomic variation in migration rates can be modeled in two alternative models:</h3>'),
				HTML('<h3><b>A. Beta(&alpha;, &beta;) distribution</b> of <i>N.m</i>:</h3>')
			)
		),

		fluidRow(
			# Sidebar panel for inputs
			column(width=4,
				# Input: Slider for the number of bins
				sliderInput(inputId = "M", label = h3("Effective migration rates (number of migrants per generation):"), min = 0, max = 10, value = 5, step=0.1),
				sliderInput(inputId = "alphaM", label = HTML('<h3>Shape parameter &alpha;:</h3>'), min = 1, max = 20, value = 10, step=0.1),
				sliderInput(inputId = "betaM", label = HTML('<h3>Shape parameter &beta;:</h3>'), min = 1, max = 20, value = 3, step=0.1)	
			),
			
			# Main panel for displaying outputs
			column(width=8,
				# Output: density
				plotOutput(outputId = "migration_hetero_beta")
			)
		),

		fluidRow(
			column(
				hr(),
				width = 12,
				HTML('<h3><b>B. bimodal distribution</b> of <i>N.m</i>:</h3>')
			)
		),

		fluidRow(
			# Sidebar panel for inputs
			column(width=4,
				# Input: Slider for the number of bins
				sliderInput(inputId = "Mexample", label = h3("Effective migration rates (number of migrants per generation):"), min = 0, max = 10, value = 5, step=0.1),
				sliderInput(inputId = "nLociExample", label = h3("Number of studied loci:"), min = 20, max = 1000, value = 100, step=1),
				sliderInput(inputId = "propBarrierExample", label = h3("Proportion of barriers (in %):"), min = 0, max = 100, value = 50, step=1)
			),
			
			# Main panel for displaying outputs
			column(width=8,
				# Output: density
				plotOutput(outputId = "migration_hetero_bimodal")
			)
		)

	),
	
#	box(title = h2("Model comparisons for 2 populations/species"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
	boxPlus(title = h2("Model comparisons for 1 population"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
		htmlOutput("model_comparisons_1pop"),
		HTML('<h3><b>DILS</b> performs hierarchical model comparisons.</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.</b> comparison between all models with <b>Expansion</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>}) <i>versus</i> <b>Constant size</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>}) <i>versus</i> <b>Contraction</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>})</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.</b> the last step of the model comparison is to determine whether effective size (<b><i>Ne</i></b>) is homogeneously or heterogenously distributed in genomes.</h3>')
	),
	
	boxPlus(title = h2("Model comparisons for 2 populations/species"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
		htmlOutput("model_comparisons_2pop"),
		HTML('<h3><b>DILS</b> performs hierarchical model comparisons.</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.</b> comparison between all models with <b>current isolation</b> ({SI; AM} x {<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>} x {<i>M</i><sub>homo</sub>; <i>M</i><sub>hetero</sub>}) <i>versus</i> <b>ongoing migration</b> ({IM; SC} x {<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>} x {<i>M</i><sub>homo</sub>; <i>M</i><sub>hetero</sub>})</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.a. if current isolation -></b> comparison between <b>SI</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>}) <i>versus</i> <b>AM</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>} x {<i>M</i><sub>homo</sub>; <i>M</i><sub>hetero</sub>})</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.b. if ongoing migration -></b> comparison between <b>IM</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>} x {<i>M</i><sub>homo</sub>; <i>M</i><sub>hetero</sub>}) <i>versus</i> <b>SC</b> ({<i>Ne</i><sub>homo</sub>; <i>Ne</i><sub>hetero</sub>} x {<i>M</i><sub>homo</sub>; <i>M</i><sub>hetero</sub>})</h3>'),
		HTML('<h3><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.</b> the last step is to determine whether effective size (<b><i>Ne</i></b>; panel <b>B</b>) and migration rates (<b><i>N.m</i></b>; panel <b>C</b>) are homogeneously (<b>o</b>) or heterogenously (<b>e</b>) distributed in genomes.</h3>')
	),

	boxPlus(title = HTML('<h2>Architecture of DILS</h2>'), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
		column(width=12,
			h3("FastABC is composed of two elements:"),
			h3(strong("1."), "A web interface developed in", a(span(strong("Shiny,"), style = "color:teal"), href="https://shiny.rstudio.com/", target="_blank"), "which will execute..."),
			h3(strong("2."), "...a workflow managed by", a(span(strong("Snakemake."), style = "color:teal"), href="https://snakemake.readthedocs.io/en/stable/", target="_blank")),
			htmlOutput("welcome_picture"),
			hr(),
			h3("The code is fully open-source, freely distributed on", a(span(strong("GitHub"), style = "color:teal"), href="https://github.com/popgenomics/DILS", target="_blank"), "and can be immediately redeployed on any cluster using", a(span(strong("SLURM"), style = "color:teal"), href="https://slurm.schedmd.com/documentation.html", target="_blank"), "thanks to a", a(span(strong("Singularity"), style = "color:teal"), href="https://sylabs.io/docs/#doc-3.2", target="_blank"), "image."),
			
			h3("This redeployment allows the user to modify the models to be compared, to add summary statistics, etc..."),
			h3("The workflow can be simply executed from the command line without going through the web interface."),
			
			h3("However, if desired, the web interface can also be freely hosted and linked to any cluster.")
			)
		)
	)
)

# ABC tuto
ABC_tuto <- fluidPage(
	box(title = h2("ABC"), width = NULL, status = "primary", solidHeader = TRUE,
		fluidRow(
			column(width=12,
				h3("The ABC section is organized in 5 steps:"),
				h3("4 configuration."),
				infoBox("Upload data", "pouet", icon = icon("cloud-upload"), color='navy'),
				infoBox("Data filtering", "pouet", icon = icon("bath"), color='navy'),
				infoBox("Populations/species", "pouet", icon = icon("users-cog"), color='navy'),
				infoBox("Prior distributions", "pouet", icon = icon("dice"), color='navy')
			)
		),
		hr(),
		fluidRow(
			column(width=12,
				h3("1 of execution."),
				infoBox("Run ABC", "pouet", icon = icon("microchip"), color='navy') # #556270
			)
		)
	)
)

# upload
upload_data <- fluidPage(
	tags$head(tags$style(HTML("a {color: black}"))),
	fluidRow(
		boxPlus(title = h2("Number of ABC analysis to run"), width = 6, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
		shinyjs::useShinyjs(),
		selectInput("number_of_ABC", label = h4("1 to 5 ABC analyses can be performed from the same input file"), choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5), selected = 1),
		HTML('<h4>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>DILS</b> runs freely on a computer server.<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To avoid saturating it, we have limited the number of analyses carried out at a given time and for a given input file to <b>5</b>.<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Beyond 5, you will have to upload it again whenever you want.<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The selected number of analysis cannot be modified once the choices had been checked/validated</h4>')
	),
	
	boxPlus(title = h2("Email address"), width = 6, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
		textInput("mail_address", label = h4("address to receive the download link of the results"), value = "user@gmail.com"),
		hr(),
		h4("This address will only be used for 2 things:"),
		h4(strong("1)"), "send the results of DILS to the user"),
		h4(strong("2)"), "contact users on the day when a collaborative meta-analysis will be considered")
		)
	),

	fluidRow(NULL, soldHeader = TRUE, status ="danger",
		boxPlus(title = h2("Input file upload"), height = 200,	width = 6, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
		fileInput("infile", label = NULL, accept = c('.fasta', '.fas', '.fa')),
		tags$style(".progress-bar {background-color: #1e2b37;}")
		),
	
		boxPlus(title = h2("Genomic regions"), height = 200,	width = 6, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
			prettyRadioButtons("region", label = NULL, shape = "round", status = "warning", fill = TRUE, inline = TRUE, animation = "pulse", bigger = TRUE,
			selected = "noncoding", choices = list("coding" = "coding", "non coding" = "noncoding"))
		)
	),	
	
	fluidRow(
		boxPlus(
			title = h2("Input file format"), width = 12, icon = NULL, solidHeader = TRUE, background = NULL,
			boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
			enable_label = TRUE, label_text = "CLICK TO DISPLAY AN EXAMPLE OF INPUT FILE", label_status = "success",
			h3(strong("Fasta file")),
			h3("A single fasta file containing all sequences obtained from all populations/species, and for all genes is the only inputfile to upload."),
			h3("Even sequences obtained from non-studied species can be included in the file. The user will specify the names of the species to consider after the upload ."),
			h3("Its format is largely inspired by the output of", a(strong("Reads2snp"), href="https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools", target="_blank"), " but can be post-produced without Reads2snp."),
			hr(),
			h3("All sequences's have to respect the following structure:"),
			h3(strong(">gene|species or population|individual|allele1 or allele2")),
			h3(strong("GTGATGCGTGTAGTCATG")),
			h3("With missing data only encoded by 'N'"),
			br(),
			h3(strong("Example:")),
			p(">Hmel210004_196|chi|chi.CJ560|allele1"),
			p("NNNNNNNGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|chi|chi.CJ560|allele2"),
			p("NNNNNNNGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|chi|chi.CJ564|allele1"),
			p("NTGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|chi|chi.CJ564|allele2"),
			p("NTGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATATTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|flo|flo.CS2338|allele1"),
			p("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
			p(">Hmel210004_196|flo|flo.CS2338|allele2"),
			p("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
			p(">Hmel210004_196|flo|flo.CS2341|allele1"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|ros|ros.CJ2071|allele1"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|ros|ros.CJ2071|allele2"),
			p("ATGTCTCGGCCAGTGTTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCTGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|ros|ros.CJ531|allele1"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCTGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|num|nu_sil.MJ09-4125|allele1"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACCTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTGGAATTATACAAA"),
			p(">Hmel210004_196|num|nu_sil.MJ09-4125|allele2"),
			p("ATGTCTCGGCCAGTATTATCTACACACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTGGAATTATACAAA"),
			p(">Hmel210004_196|num|nu_sil.MJ09-4184|allele1"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAGTGGAATTTTCGTCGAATTATACAAA"),
			p(">Hmel210004_196|num|nu_sil.MJ09-4184|allele2"),
			p("ATGTCTCGGCCAGTATTATCTACGCACGTGTTAGACACTTCNACTGGTCAGCCAGGAAATGGAATTTTCGTGGAATTATACAAA"),
			p(">Hmel219015_26|chi|chi.CAM25091|allele1"),
			p("GGAAATNNAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTNTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACANGAGAAAATAA"),
			p(">Hmel219015_26|chi|chi.CAM25091|allele2"),
			p("GGAAATNNAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTNTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTAAAATTTCACANGAGAAAATAA"),
			p(">Hmel219015_26|chi|chi.CAM25137|allele1"),
			p("GGAAATNNAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCANAGGAGAAAATAA"),
			p(">Hmel219015_26|chi|chi.CAM25137|allele2"),
			p("GGAAATNNAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTAAAATTTCANAAGAGAAAATAA"),
			p(">Hmel219015_26|flo|flo.CS12|allele1"),
			p("GGAAATGAAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAAGAGAAAATAA"),
			p(">Hmel219015_26|flo|flo.CS12|allele2"),
			p("GGAAATGAAAACTTTTGTATCAAGTTTGTTACGGCGATTTCGCCTAGAAGCTGGAACAAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAAGAGAAAATAA"),
			p(">Hmel219015_26|flo|flo.CS13|allele1"),
			p("NNNNNNGAAAACTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
			p(">Hmel219015_26|ros|ros.CAM1841|allele1"),
			p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
			p(">Hmel219015_26|ros|ros.CAM1841|allele2"),
			p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
			p(">Hmel219015_26|ros|ros.CAM1880|allele1"),
			p("NNNNNNNNNNNNNNNNNTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
			p(">Hmel219015_26|ros|ros.CAM1841|allele1"),
			p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
			p(">Hmel219015_26|ros|ros.CAM1841|allele2"),
			p("NNNNNNNNNNNCTTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACNNNNNNNNNNNN"),
			p(">Hmel219015_26|ros|ros.CAM1880|allele1"),
			p("NNNNNNNNNNNNNNNNNTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTTAAATTTCACAGGAGAAAATAA"),
			p(">Hmel219015_26|num|nu_sil.MJ09-4125|allele1"),
			p("NNNNNNGNNNNCTTTTGTATNAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGATCTGGTCTTCCGCACTGATATTATATTGCGAACTATGGGACAACCAATTTACGTNAAATTTCACAAGAGAAAATAA"),
			p(">Hmel219015_26|num|nu_sil.MJ09-4125|allele2"),
			p("NNNNNNTNNNNCTTTTGTATNAATTCTGTTGAGGCGATTTCGTCTAGAAGCTGTAACGAAGCCATCTGACCTGGTCTTTCGCACTGATATTATATTACGAACTATTGGACAACCAGTGTACGTNAAATTTCACAAAAGAAAATAA"),
			p(">Hmel219015_26|num|nu_sil.MJ09-4184|allele1"),
			p("NGANATGAAAACNTTTGTATCAAGTGTGTTACGGCGATTTCGTCTAGAAGCTGTAACNAAGCCATCTGATCTGGTCTTCCGCACTGATATTNTATTGCGAACTATGGGACAACCAATTTACGTNAAATTTCACANNAGAAAATAA"),
			p(">Hmel219015_26|num|nu_sil.MJ09-4184|allele2"),
			p("NGANATTAAAACNTTTGTATCAATTCTGTTGAGGCGATTTCGTCTAGAAGCTGTAACNAAGCCATCTGATCTGGTCTTTCGCACTGATATTNTATTACGAACTATTGGACAACCAGTGTACGTNAAATTTCACANNAGAAAATAA"),
			p("etc ..."),
			h3("Two genes are displayed in this example, they are named: ", strong("Hmel210004_196"), " and", strong("Hmel219015_26.")),
			h3("Four populations are present in this example, named: ", strong("chi, flo, ros and num.")),
			h3("Only species whose names are specified in the ", strong("Populations/species"), " menu are considered, but the uploaded file can contain other species."),
			h3("Two diploid individuals are sequenced for each species/population. For example for chi: chi.CJ560 and chi.CJ564. This number can obviously vary between species/populations, according to the sequencing strategy and its success.")
		)
	),
	
	fluidRow(align="left",
		boxPlus(title = h2("Information extracted from the uploaded file"), width = 6, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
			uiOutput("upload")
		),
		
		boxPlus("", width = 6, solidHeader = TRUE, status = "info",
			prettyCheckbox(inputId = "check_upload", shape = "round", value = FALSE,
			label = strong("Please check/valid your choices"), icon = icon("check"),
			animation = "tada", status = "success", bigger = TRUE)
		)
	)


)


filtering <- fluidPage(
	fluidRow(
		column(width = 4,
			boxPlus(title = h2("Maximum proportion of missing data (N, gaps, ...)"), height = 225, width = NULL, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				sliderInput("max_N_tolerated", label = NULL,	min = 0, max = 1, value = 0.1, step = 0.005)
			),
			boxPlus(
				title = h3("max_N_tolerated"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "INFORMATION", label_status = "primary",
				h3("-Float between 0.0 and 1.0."),
				h3("-Defines the maximum proportion of N in the sequence of a gene beyond which this sequence is not considered."),
				hr(),
				h3(a(span(strong("Example", style = "color:blue")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/max_N_tolerated.png", target="_blank"))
			)
		),
		
		column(width = 4,
			boxPlus(title = h2("Minimum sequence length per gene"), height = 225, width = NULL, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				numericInput("Lmin", label = NULL, value = 30, min = 1, max = 10000)
			),

			boxPlus(
				title = h3("Lmin"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "INFORMATION", label_status = "success",
				
				h3("-Positive integer (>0)"),
				h3("-Minimum number of treatable sites below which a gene is removed from the analysis."),
				br(),
				h3("In a noncoding sequence: a site is an alignment of nucleotides for a single given nucleotide position, including all individuals among the species considered."),
				br(),
				h3("In a coding sequence: a site is an alignment for a given codon, comprising all the individuals among the species considered."),
				h3("A coding position is not considered if:"),
				h3("	-a codon alignment contains a non-synonymous polymorphism."),
				h3("	-more than two codons segregate (even synonyms)."),
				h3("	-at least one N is found in a codon, in an individual."),
				br(),
				h3("Number of positions to consider = (number of ", strong("monomorphic positions"), "that can be considered) + (number of ", strong("biallelic positions"), "that can be considered)."),
				hr(),
				h3(a(span(strong("Example", style = "color:green")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/Lmin.png", target="_blank"))
			)
		),
		
		column(width = 4,
			boxPlus(title = h2("Minimum number of sequences per gene and per population/species"), height = 225, width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				numericInput("nMin", label = NULL, value = 12, min = 2)
			),
			
			boxPlus(
				title = h3("nMin"), width = NULL, icon = NULL, solidHeader = TRUE, background = NULL,
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "INFORMATION", label_status = "warning",
				h3(strong("DILS"), " starts for each gene by eliminating individual sequences containing too many N and gaps", span(strong("(max_N_tolerated; blue box)", style = "color:blue")), "."),
				br(),
				h3("If for a gene and", strong("within a population/species"), "there are fewer than", strong("nMin"), "sequences left, then the gene is not considered in the ", strong("ABC"), "analysis."),
				br(),
				h3(strong("If an outgroup is specified:")),
				h3(strong("nMin"), " becomes the number of sequences sampled for each species at each locus, to produce a standardized joint site frequency spectrum (jSFS) used by ", strong("ABC"), "."),
				hr(),
				h3(a(span(strong("Example", style = "color:orange")), href="https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/nMin.png", target="_blank"))
			)
		)
	),
	
	fluidRow(
		box("", width = 12, solidHeader = TRUE, status = "info",
			prettyCheckbox(inputId = "check_filtering", shape = "round", value = FALSE,
				label = strong("Please check/valid your choices"), icon = icon("check"),
				animation = "tada", status = "success", bigger = TRUE)
		)
	)
)


populations <- fluidPage(
	fluidRow(
		boxPlus(title = h2("Number of populations/species"), height = NULL, width = 4, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
		#	### FOR THE MOMENT : ONLY ONE POSSIBLE CHOICE --> 2 POPULATIONS SPECIES.
		#	### UNCOMMENT THE NEXT TWO LINES WHEN THE ANALYSES FOR 1 OR 4 POPULATIONS WILL BE READY
			#	prettyRadioButtons("nspecies", label = h3("Number of gene pools"), shape = "round", status = "primary", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
			#	choices = list("One gene pool" = 1, "Two gene pools" = 2, "Four gene pools" = 4), selected = 2),

			prettyRadioButtons("nspecies", label = NULL, shape = "round", status = "primary", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
			choices = list("One gene pool" = 1, "Two gene pools" = 2), selected = 2),
			HTML('<h4>If <b>set to one</b>: a model for a single panmictic population is evaluated.<br>If <b>set to two</b>: a model of divergence between two populations/species is evaluated.<br></h4>'),
			em(strong(h4('Analysis for 4 populations/species will be soon available'))),
			uiOutput("input_names_ui")
		),
		
		boxPlus(title = h2("Presence of an outgroup"), height = NULL, width = 4, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
			prettyRadioButtons("presence_outgroup", label = NULL, shape = "round", status = "warning", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
			choices = list("no" = "no", "yes" = "yes"), selected = "no"),
			uiOutput("input_names_outgroup_ui"),
			HTML('<h4>If <b>set to yes</b>: the used joint site frequency spectrum (jSFS) is <b>unfolded</b>, and the mutation rate for each locus is corrected by its divergence with the outgroup.</h4>'),
			HTML('<h4>If <b>set to no</b>: the used jSFS is <b>folded</b>, and the mutation rate is the same for all loci.</h4>')
		),
		uiOutput("size_change")
	),

	fluidRow(
		box("", width = 12, solidHeader = TRUE, status = "info",
			prettyCheckbox(inputId = "check_populations", shape = "round", value = FALSE,
			label = strong("Please check/valid your choices"), icon = icon("check"),
			animation = "tada", status = "success", bigger = TRUE)
		)
	)
)


prior <- fluidPage(
	fluidRow(
		column(width = 6,
			uiOutput("prior_mutation"),
		
			uiOutput("prior_Ne")
		),
		
		column(width = 6,
			uiOutput("prior_times"),
		
			uiOutput("prior_migration")
		)
	),
		
	fluidRow(
		box("", width = 12, solidHeader = TRUE, status = "info",
			prettyCheckbox(inputId = "check_prior", shape = "round", value = FALSE,
			label = strong("Please check/valid your choices"), icon = icon("check"),
			animation = "tada", status = "success", bigger = TRUE)
		)
	)
)



run_ABC <- fluidPage(
	# PRINT INFOX BOXES OF CHECKING
	fluidRow(
		boxPlus(
			width = 12,
			uiOutput('check_upload_info'),
			uiOutput('check_filtering_info'),
			uiOutput('check_populations_info'),
			uiOutput('check_prior_info')
		)
	),
	
	# PRINT INPUT
	fluidRow(
		column(width = 12,
			boxPlus(
				title = h2("Information summary"), width = NULL, icon = "fa fa-heart", solidHeader = TRUE, gradientColor = "teal",
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "Please check the following information", label_status = "success",

				tableOutput("parameters")
			)
		),
		
		# RUN ABC
		uiOutput("run_ABC")
	)
)


upload_results <- fluidPage(
	# upload results
	fluidRow(NULL, soldHeader = TRUE, status ="danger",
		boxPlus(title = h2("Results to upload (i.e, DILS's archived output)"), height = 240,	width = 12, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
			HTML('<h3>This is the DILS&#39;s archived output that you downloaded from the link sent to your email address</h3>'),
			fileInput("results", label = NULL, accept=c('.tar.gz')),
			tags$style(".progress-bar {background-color: #1e2b37;}")
		)
	),
	uiOutput("observed_columns_to_display"),
	uiOutput("display_uploaded_results")
)


user_dataset <- fluidPage(
	uiOutput("visualization_data")
)


collaborative <- fluidPage(
	tags$head(tags$script('
		var dimension = [0, 0];
		$(document).on("shiny:connected", function(e) {
		dimension[0] = window.innerWidth;
		dimension[1] = window.innerHeight;
		Shiny.onInputChange("dimension", dimension);
		});
		$(window).resize(function(e) {
		dimension[0] = window.innerWidth;
		dimension[1] = window.innerHeight;
		Shiny.onInputChange("dimension", dimension);
		});
	')),
	
	h3(''),
	hr(),
	
	uiOutput("page_greyzone")
)


information <- fluidPage(
	fluidRow(
		column(width = 12,
			#box(title = h2("Citations"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
			boxPlus(title = h2("Citations"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				h3("Please, in case of publication of a study using DILS, do not forget to quote the following references:"),
				h3(code('Csilléry, Katalin, Olivier François, and Michael GB Blum. "abc: an R package for approximate Bayesian computation (ABC)." Methods in ecology and evolution 3.3 (2012): 475-479.')),
				h3(code('Pudlo, Pierre, Jean-Michel Marin, Arnaud Estoup, Jean-Marie Cornuet, Mathieu Gautier, and Christian P. Robert. "Reliable ABC model choice via random forests." Bioinformatics 32, no. 6 (2015): 859-866.')),
				h3(code('Roux, Camille, Christelle Fraisse, Jonathan Romiguier, Yoann Anciaux, Nicolas Galtier, and Nicolas Bierne. "Shedding light on the grey zone of speciation along a continuum of genomic divergence." PLoS biology 14, no. 12 (2016): e2000234.'))
			),
			
			#box(title = h2("Acknowledgment"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
			boxPlus(title = h2("Acknowledgment"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				h3("Please, if you use this online version of DILS, do not forget to recognize and acknowledge the free provision of calculation cores by France Bioinformatique"),
				h3(code("The demographic inferences were conducted on the IFB Core Cluster which is part of the National Network of Compute Resources (NNCR) of the", a(span(strong("Institut Français de Bioinformatique (IFB)."), style = "color:teal"), href="https://www.france-bioinformatique.fr/fr", target="_blank")))
			),
			
			#box(title = h2("Partners"), width = 12, solidHeader = TRUE, background = NULL, status = "primary",
			boxPlus(title = h2("Partners"), width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				mainPanel(htmlOutput("logos"))
			)
		)
	)
)

ui <- dashboardPage(
	
	#skin = "black",
	dashboardHeader(title = "menu DILS",
#		tags$li(class = "dropdown", socialButton(url = "https://github.com/popgenomics/ABConline", type = "github"), tags$img(height = "auto"))
		tags$li(class="dropdown", tags$a(href="https://github.com/popgenomics/DILS", icon("github"), "Source Code", target="_blank")),
		tags$li(class="dropdown", tags$a(href="https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection", icon("envelope"), "Help/Discussion", target="_blank"))

	),

	dashboardSidebar(
		sidebarMenu(
#			style = "position: fixed; overflow: visible; width: inherit",
			menuItem(("Welcome"), tabName = "welcome", icon = icon("door-open", class="door-open")),

			menuItem(('ABC'), tabName = "ABC", icon = icon("industry"),
				menuSubItem(("Upload data"), tabName = "upload", icon = icon("cloud-upload")),
				menuSubItem(("Data filtering"), tabName = "filtering", icon = icon("bath")),
				menuSubItem(("Populations/species"), tabName = "populations", icon = icon("users-cog")),
				menuSubItem(("Prior distributions"), tabName = "simulations", icon = icon("dice")),
				menuSubItem(("Run ABC"), tabName = "run_abc", icon = icon("microchip"))),

			menuItem(("Results visualization"), tabName = "Results_visualization", icon = icon("chart-pie"),
				menuSubItem(("Upload results"), tabName = "upload_results", icon = icon("cloud-upload")),
				menuSubItem(("User's dataset"), tabName = "user_dataset", icon = icon("lock")),
				menuSubItem(("Collaborative science"), tabName = "collaborative", icon = icon("lock-open"))
			),
			
			menuItem(("Information"), tabName = "information", icon = icon("info-circle"))
		)
	),
	
	dashboardBody(
		# tags
		## bar du menu sur le cote
		## items du menu quand on passe dessus a la souris
		## background du body
		## background du text du menu quand on passe dssus à la souris
		tags$head(tags$style(HTML('
			/* HEADER */
			/* fond derriere le nom du header;  nom du header */
			.skin-blue .main-header .logo { background-color: #1e2b37; color: #ffffff; font-size: 24px; height: 50px;}

			/* couleur de fond du header sous la souris */
			.skin-blue .main-header .logo:hover { background-color: #1e2b37; height: 50px;}

			/* toute la partie droite de la barre du header */
			.skin-blue .main-header .navbar { background-color: #1e2b37; height: 40px; font-size: 20px}

			/* toute la partie droite de la barre du header sous la souris*/
			.skin-blue .main-header .navbar:hover{ background-color: #1e2b37; font-color: #C7F464; height: 40px; font-size: 20px}
			
			/* bouton menu dans le header: background et petits traits	*/
			.skin-blue .main-header .navbar .sidebar-toggle{ background-color: #1e2b37; color: #ffffff; height: 40px; }
	
			/* bouton menu sous la souris dans le header: background et petits traits */	
			.skin-blue .main-header .navbar .sidebar-toggle:hover{ background-color: #1e2b37 ;color: #C7F464; height: 40px; }
			
			/* SIDEBAR */
			/* taille de la police */
			.main-sidebar { font-size: 18px; }
			
			/* couleur du fond du menu */
			.skin-blue .main-sidebar { background-color: #556270;font-size: 18px; }
			
			/* couleur des elements du menu sous la souris */
			.skin-blue .sidebar-menu>li.active>a, .skin-blue .sidebar-menu>li:hover>a { background-color: #1e2b37;font-size: 18px; }
			
			/* couleur du texte du menu quand on passe la souris dessus */	
			.skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{ color: #C7F464;font-size: 18px; }
			
			/* elements du menu selectionnes */
			.skin-blue .main-sidebar .sidebar .sidebar-menu .active a{ background-color: #1e2b37; color: #C7F464;font-size: 18px; }
		
			/* other links in the sidebarmenu */
			.skin-blue .main-sidebar .sidebar .sidebar-menu a{color: #ffffff;font-size: 18px;}
	
			/* BODY */
			/* couleur du background du body */
			.content-wrapper, .right-side { background-color: ghost-white; }
			
			/* TABSET */
			.tabbable > .nav > li > a {background-color: #ffffff; color:#556270;font-size: 18px;}
			.tabbable > .nav > li[class=active] > a {background-color: #556270; color:#C7F464;font-size: 18px;}
			
			/* FILEINPUT */
			.btn-file { background-color:#556270; border-color: color:#C7F464; color:#C7F464; font-size: 18px; }
			/*.btn-file { background-color:#556270; border-color: color:#C7F464; color:#C7F464;}*/

			/* BUTTONS */
/*			.btn .btn-social-icon .btn-github a{color:#fff; background:#556270; height: 100px; font-size: 30px;}*/
			
			/* BOX */
/*			.box.box-solid.box-primary>.box-header { color:#fff; background:#556270; font-size: 0px; } */
/*			.box.box-solid.box-primary{border-bottom-color:#556270; border-left-color:#556270; border-right-color:#556270; border-top-color:#556270;} */

		'))),
		
		setShadow(class = "box"),
#		shinyDashboardThemes(
			#theme = "boe_website"
#			theme = "poor_mans_flatly"
#		),
	
		tabItems(
			# Welcome
			tabItem(tabName = "welcome",
				welcome_page
			),
			
			tabItem(tabName = "ABC",
				ABC_tuto
			),
	
			# Upload data
			tabItem(tabName = "upload",
				upload_data
			),
			
			#	Filtering
			tabItem(tabName = "filtering",
				filtering
			),
			
			# Populations
			tabItem(tabName = "populations",
				populations
			),
		 
			# Simulations
			tabItem(tabName = "simulations",
				prior
			),

			# Run the ABC inferences
			tabItem(tabName = "run_abc",
				run_ABC
			),
			
			# Upload the DILS's results
			tabItem(tabName = "upload_results",
				upload_results
			),
			
			# Plot the distributions of statistics
			tabItem(tabName = "user_dataset",
				user_dataset
			),
			
			# Plot the distributions of statistics
			tabItem(tabName = "collaborative",
				collaborative
			),
			
			# Information
			tabItem(tabName = "information",
				information
			)
		)
	)
)


# Define server logic required to draw a histogram
server <- function(input, output, session = session) {
	options(shiny.maxRequestSize=4000*1024^2)

	#	WELCOME
	output$welcome_picture <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/dag_2pops.pdf.png align="middle" height="auto" width="100%" margin="0 auto">'
		)
		}
	)
	
	output$overview_DILS_picture <-
		renderText({
			c('<img src=https://github.com/popgenomics/ABConline/blob/master/webinterface/pictures_folder/overview.png?raw=true align="middle" height="auto" width="50%" margin="0 auto">')
		}
	)
	
	output$models_picture_1pop <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/models_onePop.png align="middle" height="auto" width="30%" margin="0 auto">'
		)
		}
	)
	
	output$models_picture_2pop <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/models_2pops.png align="middle" height="auto" width="30%" margin="0 auto">'
		)
		}
	)
	
	output$models_picture_4pop <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/models_4pops.png align="middle" height="auto" width="50%" margin="0 auto">'
		)
		}
	)
	
	output$model_comparisons_1pop <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/figure_2.png align="middle" height="auto" width="50%" margin="0 auto">'
		)
		}
	)
	
	output$model_comparisons_2pop <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/figure_3.png align="middle" height="auto" width="50%" margin="0 auto">'
		)
		}
	)
	
	output$homo_hetero <-
		renderText({
		c(
		'<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/homo_hetero.png align="middle" height="auto" width="50%" margin="0 auto">'
		)
		}
	)
	
	## get the directory
	global <- reactiveValues(datapath = getwd()) 
	
	## example of genomic heterogeneity
	output$genomic_hetero <- renderPlot({
		par(las=1)
		y_points = dbeta(0:100/100, input$alpha, input$beta)
		x_points = 0:100/100 * input$Ne/( input$alpha / (input$alpha + input$beta) )
		plot(x_points, y_points, type='l', xlab = expression(paste("Genomic distribution of ", italic('Ne'), sep=" ")), ylab='density', main=expression(italic("Example of genomic distributions that DILS will try to infer")), col="white", cex.main = 1.5, cex.axis = 1.5, cex.lab=1.5, xlim=c(min(c(x_points, input$Ne*1.2)), max(c(x_points, input$Ne*1.2))))
		
		x_points = c(0, x_points, max(x_points), 0)
		y_points = c(0, y_points, 0, 0)
		polygon(x_points, y_points, border = 'NA', col=viridis_pal(option="D")(2)[2])
		abline(v=input$Ne, lwd=8, col=viridis_pal(option="D")(2)[1])
	})
	
	output$migration_hetero_beta <- renderPlot({
		par(las=1)
		y_points = dbeta(0:100/100, input$alphaM, input$betaM)
		x_points = 0:100/100 * input$M/( input$alphaM / (input$alphaM + input$betaM) )
		plot(x_points, y_points, type='l', xlab = expression(paste("Genomic distribution of ", italic('N.m'), sep=" ")), ylab='density', main=expression(italic("Beta(&alpha;, &beta;) distribution")), col="white", cex.main = 1.5, cex.axis = 1.5, cex.lab=1.5, xlim=c(min(c(x_points, input$M*1.2)), max(c(x_points, input$M*1.2))))
		
		x_points = c(0, x_points, max(x_points), 0)
		y_points = c(0, y_points, 0, 0)
		polygon(x_points, y_points, border = 'NA', col=viridis_pal(option="D")(2)[2])
		abline(v=input$M, lwd=8, col=viridis_pal(option="D")(2)[1])
	})
	
	output$migration_hetero_bimodal <- renderPlot({
		nLociExample = input$nLociExample
		barriers = as.integer(input$propBarrierExample/100 * nLociExample)
		distribution = c(rep(0, barriers), rep(input$Mexample, nLociExample-barriers))
		hist(distribution, xlab = expression(paste("Genomic distribution of ", italic('N.m'), sep=" ")), ylab='Number of loci', main=expression(italic("Bimodal distribution of N.m")), col=viridis_pal(option="D")(2)[2], cex.main = 1.5, cex.axis = 1.5, cex.lab=1.5, xlim=1.2*range(distribution))
	})

	#	UPLOAD DATA
	## GET THE SUMMARY STATS ABOUT THE UPLOADED FILE
	## list of species
	list_species = reactive({
		if(is.null(input$infile)){return ()}
		withProgress(message = 'Getting the species', detail = NULL, value = 0, {
		incProgress(1/2)
		print(input$infile$datapath)
		return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f2 | sort -u", sep=" "), intern = T))
		incProgress(1/2)
		})
	})
	output$list_species <- renderDataTable(data.frame("species" = list_species()))

	## list of individuals
	list_individuals = reactive({
		if(is.null(input$infile)){return ()}
		withProgress(message = 'Getting the individuals', detail = NULL, value = 0, {
		incProgress(1/2)
		return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f3 | sort -u", sep=" "), intern = T))
		incProgress(1/2)
		})
	})
	output$list_individuals <- renderDataTable(data.frame("individuals" = list_individuals()))

	## list of loci
	list_loci = reactive({
		if(is.null(input$infile)){return ()}
		withProgress(message = 'Getting the loci', detail = NULL, value = 0, {
		incProgress(1/2)
		return(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f1 | sort -u | cut -d'>' -f2", sep=" "), intern = T))
		incProgress(1/2)
		})
	})
	output$list_loci <- renderDataTable(data.frame("loci" = list_loci()))

	## Summary Stats code ##
	# this reactive output contains the summary of the dataset and display the summary in table format
	nSpecies = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f2 | sort -u", sep=" "), intern = T))})
	nIndividuals = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f3 | sort -u", sep=" "), intern = T))})
	nLoci = reactive({length(system(paste("cat", input$infile$datapath, "| grep '>' | cut -d '|' -f1 | sort -u", sep=" "), intern = T))})
	
	output$general_information <- renderTable({
		if(is.null(input$infile)){return ()}
		withProgress(message = 'Producing the table of information', detail = NULL, value = 0, {
		incProgress(1/2)
		return(data.frame("nSpecies" = length(list_species()), "nIndividuals" = length(list_individuals()), "nLoci" = length(list_loci())))
		incProgress(1/2)
		})
	})

	## MainPanel tabset renderUI code ##
	# the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
	# Until the file is loaded, app will not show the tabset.
	output$upload <- renderUI({
		if(is.null(input$infile)) {return(loadingState())}
		else
			tabsetPanel(
				tabPanel("General information", tableOutput("general_information")),
				tabPanel("List of individuals", dataTableOutput("list_individuals")),
				tabPanel("List of populations or species", dataTableOutput("list_species")),
				tabPanel("List of loci", dataTableOutput("list_loci"))
			)
	})

	## Number of ABC analysis to perform
	#	observeEvent(input$number_of_ABC_validation, {
	observe(if(input$check_upload){
		shinyjs::disable("number_of_ABC")
	})
	
	# POPULATIONS/SPECIES
	output$input_names_ui <- renderUI({
		if(is.null(input$infile)) {return()}
		nspecies = as.integer(input$nspecies)
		lapply(1:nspecies, function(i){
			selectInput(paste0("name", LETTERS[i]), label = paste0("name of species/population ", LETTERS[i]), choices = list_species(), selected = paste0("name", LETTERS[i]))
		})
	})
	
	output$input_names_outgroup_ui <- renderUI({
		if(is.null(input$infile)) {return()}
			presence_outgroup = input$presence_outgroup
		if(presence_outgroup == 'no'){
		# nameOutgroup = 'NA' # CONFIG_YAML
		}else{
			selectInput("nameOutgroup", label = "name of the outgroup species", choices = list_species())
		# nameOutgroup = 'name specified by the selectInput # CONFIG_YAML
		}
	})

	output$size_change <- renderUI({
		nspecies = as.integer(input$nspecies)
		if(nspecies==2){
		boxPlus(title = h2("Size change over time"), height = NULL, width = 4, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
			prettyRadioButtons("population_growth", label = h3("Assuming constant or variable population sizes"), shape = "round", status = "danger", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
			choices = list("constant" = "constant", "variable" = "variable"), selected = "constant"),
			HTML('<h4>If <b>set to constant</b>: the sizes of the daughter populations differ from that of the ancestral population from the split, and then they remain constant.<h4>'),
			HTML('<h4>If <b>set to variable</b>: daughter populations each have two distinct sizes, one at the origin of the population and one present since T<sub>dem</sub> generations.</h4>')
			)
		}else{
			boxPlus(title = h2("Size change over time"), height = NULL, width = 4, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				prettyRadioButtons("population_growth", label = h3("The ABC analysis will test for variation in population size"), shape = "round", status = "danger", fill = TRUE, inline = FALSE, animation = "pulse", bigger = TRUE,
				choices = list("variable" = "variable"), selected = "variable"),
				em(strong(h4('The ABC analysis will compute the probabilities of models of constant population size (a single Ne as parameter), of population expansion and of population contraction (3 parameters: current Ne, ancestral Ne, time of demographic change).'))),
				return(NULL) # return(NULL) to hide this box if nspecies==1
			)
		}
	})

	# PRINT INPUT
	output$parameters <- renderTable({
		if(is.null(input$infile)) {return()}
		mail_address = input$mail_address # email address
		config_yaml = "not to be specified" # name of the config_yaml used by snakemake
		infile = input$infile$name # name of the input fasta file
		region = input$region # coding or noncoding
		nspecies = input$nspecies # number of species to simulate
		if(nspecies == 1){
			species_names = c(input$nameA)
			species_names_row = c('nameA')
		}else{
			M_min = input$M_min
			M_max = input$M_max
			if(nspecies == 2){
				species_names = c(input$nameA, input$nameB)
				species_names_row = c('nameA', "nameB")
			} else{
				if(nspecies == 4){
					species_names = c(input$nameA, input$nameB, input$nameC, input$nameD) # name of the simulated species
					species_names_row = c('nameA', "nameB", "nameC", "nameD")
				}else{
					species_names = "NA"
					species_names_row = "NA"
				}
			}
		}
		
		if(input$presence_outgroup == 'yes'){
			nameOutgroup = input$nameOutgroup
		}else{
			nameOutgroup = "NA"
		}
		
		Lmin = input$Lmin
		nMin = input$nMin
		mu = input$mu
		rho_over_theta = input$rho_over_theta
		N_min = input$N_min
		N_max = input$N_max
		Tsplit_min = input$Tsplit_min
		Tsplit_max = input$Tsplit_max
		population_growth = input$population_growth
	
		if(nspecies == 1){
			res = matrix(c(mail_address, config_yaml, infile, region, nspecies, species_names, nameOutgroup, Lmin, nMin, mu, rho_over_theta, N_min, N_max, Tsplit_min, Tsplit_max, population_growth), ncol = 1)
			row_names_res = c("user's email address", "config_yaml", "infile", "region", "nspecies", species_names_row, "nameOutgroup", "Lmin", "nMin", "mu", "rho_over_theta", "N_min", "N_max", "Tchanges_min", "Tchanges_max", "population_growth")
		}else{
			res = matrix(c(mail_address, config_yaml, infile, region, nspecies, species_names, nameOutgroup, Lmin, nMin, mu, rho_over_theta, N_min, N_max, Tsplit_min, Tsplit_max, M_min, M_max, population_growth), ncol = 1)
			row_names_res = c("user's email address", "config_yaml", "infile", "region", "nspecies", species_names_row, "nameOutgroup", "Lmin", "nMin", "mu", "rho_over_theta", "N_min", "N_max", "Tsplit_min", "Tsplit_max", "M_min", "M_max", "population_growth")
		}
		
		row.names(res) = row_names_res
		colnames(res) = c("entries")
		
		res	
	}, rownames = TRUE, colnames = TRUE)
	
	# RUN ABC
	## only show the action button RUN ABC if a file is uploaded and 4 checkings were made
	output$run_ABC <- renderUI({
		if(is.null(input$infile)==FALSE && input$check_upload == TRUE && input$check_filtering == TRUE && input$check_populations == TRUE && input$check_prior == TRUE){
			a <-column(width = 12,
				boxPlus(
					title = h2("Run ABC"), width = NULL, icon = "fa fa-heart", solidHeader = TRUE, gradientColor = "teal",
					boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
					enable_label = TRUE, label_text = "Are you ready?", label_status = "danger",
					h3("Submission of the ABC workflow"),
					actionButton("runABC", label = "Run the ABC", size = 'md', width = '100%', fullwidth = TRUE),
					h3("Number of submitted analysis (you can change the studied populations/species and the priors between 2 analysis):"),
					verbatimTextOutput('nClicks'),
					hr(),
					h3("Timestamp of the last submitted analysis:"),
					verbatimTextOutput('time_stamp'),
					hr(),
					h3("DILS command line:"),
					verbatimTextOutput("DILS_command")
				)
			)
		}else{return()}
	})
	
	## print the number of clicks on the "run the ABC" button
	output$nClicks <- renderText({ input$runABC })
	
	## removing the "Run the ABC" button after clicking on it
	observeEvent(input$runABC, {if (input$runABC == input$number_of_ABC)	removeUI(selector='#run_ABC', immediate=TRUE)}, autoDestroy=TRUE)
	
	## get a time stamp when clicking on the "Run the ABC" button
	time_stamp <- reactiveVal(0)
	observeEvent( input$runABC, {time_stamp(system('echo $(mktemp -d -t XXXXXXXXXX | cut -d"/" -f3)', intern=T))})
	output$time_stamp <- renderText({time_stamp()})

	## write the yaml file
	output$global <- renderText({global$datapath})

	observeEvent(input$runABC, {	
		if(input$presence_outgroup == 'yes'){
			nameOutgroup = input$nameOutgroup
		}else{
			nameOutgroup = "NA"
		}
		
		yaml_name = paste(global$datapath, '/', time_stamp(), '.yaml', sep='')
		write(paste("mail_address:", input$mail_address, sep=' '), file = yaml_name, append=F)
		write(paste("infile:", input$infile$datapath, sep=' '), file = yaml_name, append=T) # TO CHECK
		write(paste("region:", input$region, sep=' '), file = yaml_name, append=T)
		write(paste("nspecies:", input$nspecies, sep=' '), file = yaml_name, append=T)
		write(paste("nameA:", input$nameA, sep=' '), file = yaml_name, append=T)
		if(input$nspecies == 2){
			write(paste("nameB:", input$nameB, sep=' '), file = yaml_name, append=T)
		}
		
		write(paste("nameOutgroup:", nameOutgroup, sep=' '), file = yaml_name, append=T)
		write(paste("config_yaml:", yaml_name, sep=' '), file = yaml_name, append=T)
		write(paste("timeStamp:", time_stamp(), sep=' '), file = yaml_name, append=T)
		if(input$nspecies == 2){
			write(paste("population_growth:", input$population_growth, sep=' '), file = yaml_name, append=T)
			write(paste("modeBarrier:", input$modeBarrier, sep=' '), file = yaml_name, append=T)
		}
		
		write(paste("max_N_tolerated:", input$max_N_tolerated, sep=' '), file = yaml_name, append=T)
		write(paste("Lmin:", input$Lmin, sep=' '), file = yaml_name, append=T)
		write(paste("nMin:", input$nMin, sep=' '), file = yaml_name, append=T)
		write(paste("mu:", input$mu, sep=' '), file = yaml_name, append=T)
		write(paste("rho_over_theta:", input$rho_over_theta, sep=' '), file = yaml_name, append=T)
		write(paste("N_min:", input$N_min, sep=' '), file = yaml_name, append=T)
		write(paste("N_max:", input$N_max, sep=' '), file = yaml_name, append=T)
		if(input$nspecies == 1){
			write(paste("Tchanges_min:", input$Tsplit_min, sep=' '), file = yaml_name, append=T)
			write(paste("Tchanges_max:", input$Tsplit_max, sep=' '), file = yaml_name, append=T)
		}
		
		if(input$nspecies == 2){
			write(paste("Tsplit_min:", input$Tsplit_min, sep=' '), file = yaml_name, append=T)
			write(paste("Tsplit_max:", input$Tsplit_max, sep=' '), file = yaml_name, append=T)
			write(paste("M_min:", input$M_min, sep=' '), file = yaml_name, append=T)
			write(paste("M_max:", input$M_max, sep=' '), file = yaml_name, append=T)
		}
	})
		
#	write.table(parameters(), yaml_name, col.names=F, row.names=F, quote=F)
	
	## snakemake command
	DILS_command <- reactiveVal(0)
	#observeEvent( input$runABC, {DILS_command(paste('snakemake -p -j 999 --snakefile ../2pops/Snakefile --configfile config_', time_stamp(), '.yaml --cluster-config ../cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}"', sep=''))})
	observeEvent( input$runABC, {DILS_command(paste('DILS_', input$nspecies, 'pop.sh ', time_stamp(), '.yaml &', sep=''))})
	output$DILS_command <- renderText({DILS_command()})
	
	## Check upload
	output$check_upload_info <- renderUI({
		if(input$check_upload == FALSE) {
			a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "red", fill = TRUE, width = 3)
		} else if(input$check_upload == TRUE){
			a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("cloud-upload"), color = "green", fill = TRUE, width = 3)
		}
	})
	
	## Check filtering
	output$check_filtering_info <- renderUI({
		if(input$check_filtering == FALSE) {
			a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("bath"), color = "red", fill = TRUE, width = 3)
		} else if(input$check_filtering == TRUE){
			a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("bath"), color = "green", fill = TRUE, width = 3)
		}
	})
	
	## Check populations
	output$check_populations_info <- renderUI({
		if(input$check_populations == FALSE) {
			a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "red", fill = TRUE, width = 3)
		} else if(input$check_populations == TRUE){
			a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("users-cog"), color = "green", fill = TRUE, width = 3)
		}
	})
	
	## Check prior
	output$check_prior_info <- renderUI({
		if(input$check_prior == FALSE) {
			a <- infoBox(title= NULL, value = h4("NON CHECKED"), subtitle = NULL, icon = icon("dice"), color = "red", fill = TRUE, width = 3)
		} else if(input$check_prior == TRUE){
			a <- infoBox(title= NULL, value = h4("CHECKED"), subtitle = NULL, icon = icon("dice"), color = "green", fill = TRUE, width = 3)
		}
	})
	
	
	#tag$style(type = 'text/css', '.tab-panel{ background-color: red; color: white}')
	## RESULT VISUALIZATION
	### user information
	# all data contained in the archive
	allData <- reactive({
		
		fileName = input$results
		if(is.null(fileName)){
			return (NULL)
		}else{
			allData = list()
			
			untar(fileName$datapath, exdir = getwd())
			rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
		
			users_infos = read.table(paste(rootName, "/general_infos.txt", sep=''), h=F, sep=',')
			hierarchical = read.table(paste(rootName, "/modelComp/hierarchical_models.txt", sep=''), h=F, sep='\t')
			ABCstatGlobal = read.table(paste(rootName, "/ABCstat_global.txt", sep=''), h=T)
			priorfile = read.table(paste(rootName, '/best_model/priorfile.txt', sep=''), h=T)
			posterior = read.table(paste(rootName, '/best_model/posterior_bestModel.txt', sep=''), h=T)
			distribution_PCA = read.table(paste(rootName, '/distribution_PCA.txt', sep=''), sep='\t', h=T)
			contribution_PCA = read.table(paste(rootName, "/table_contrib_PCA_SS.txt", sep=''), h=T, sep='\t')
			coord_PCA_SS = read.table(paste(rootName, "/table_coord_PCA_SS.txt", sep=''), h=T, sep='\t')
			eigen = read.table(paste(rootName, "/table_eigenvalues_PCA_SS.txt", sep=''), h=T, sep='\t')
			gof_table = read.table(paste(rootName, "/gof/goodness_of_fit_test.txt", sep=''), h=T)
			gof2_table = read.table(paste(rootName, "/gof_2/goodness_of_fit_test.txt", sep=''), h=T)
			gof_sfs = read.table(paste(rootName, "/gof/gof_sfs.txt", sep=''), h=T)
			gof2_sfs = read.table(paste(rootName, "/gof_2/gof_sfs.txt", sep=''), h=T)
			
			meta = read.table('metaanalysis.txt', sep='\t', h=T)

			# if 2 species
			if(users_infos[1,2]==2){
				yaml = read_yaml(paste(rootName, '/config.yaml', sep=''))
				allData[['yaml']] = yaml
				
				Nref = as.numeric(read.table(paste(rootName, '/Nref.txt', sep='')))
				allData[['Nref']] = Nref 
				
				optimized_posterior = read.table(paste(rootName, '/best_model_5/posterior_bestModel.txt', sep=''), h=T)
				allData[['optimized_posterior']] = optimized_posterior 
				
				posterior_RF = read.table(paste(rootName, '/best_model/posterior_summary_RandomForest_bestModel.txt', sep=''), h=T)
				allData[['posterior_RF']] = posterior_RF
				
				optimized_posterior_RF = read.table(paste(rootName, '/best_model_5/posterior_summary_RandomForest_bestModel.txt', sep=''), h=T)
				allData[['optimized_posterior_RF']] = optimized_posterior_RF
				
				locus_spe = read.table(paste(rootName, "/locus_modelComp/locus_specific_modelComp.txt", sep=''), h=T)
				allData[['locus_spe']] = locus_spe
				
				locus_infos = read.table(paste(rootName, "/", users_infos[2,2], "_", users_infos[3,2], "_infos.txt", sep=''), h=T)
				allData[['locus_infos']] = locus_infos
			}else{
				if(users_infos[1,2]==1){
					table_posterior = read.table(paste(rootName, '/best_model/report_', users_infos[2,2], '.txt', sep=''), h=T, skip=4)
					allData[['table_posterior']] = table_posterior
					
					table_optimized_posterior = read.table(paste(rootName, '/best_model_7/report_', users_infos[2,2], '.txt', sep=''), h=T, skip=4)
					allData[['table_optimized_posterior']] = table_optimized_posterior
					
					optimized_posterior = read.table(paste(rootName, '/best_model_7/posterior_bestModel.txt', sep=''), h=T)
					allData[['optimized_posterior']] = optimized_posterior
					
					locus_spe = read.table(paste(rootName, "/ABCstat_loci.txt", sep=''), h=T)
					allData[['locus_spe']] = locus_spe
					
					locus_infos = read.table(paste(rootName, "/", users_infos[2,2], "_infos.txt", sep=''), h=T)
					allData[['locus_infos']] = locus_infos
				}
			}
		
			system(paste('rm -rf ', rootName, sep=''))
			
			allData[['users_infos']] = users_infos
			allData[['hierarchical']] = hierarchical
			allData[['ABCstatGlobal']] = ABCstatGlobal 
			allData[['priorfile']] = priorfile 
			allData[['posterior']] = posterior 
			allData[['list_parameters']] = colnames(posterior) 
			allData[['distribution_PCA']] = distribution_PCA 
			allData[['contribution_PCA']] = contribution_PCA 
			allData[['coord_PCA_SS']] = coord_PCA_SS
			allData[['eigen']] = eigen 
			allData[['gof_table']] = gof_table 
			allData[['gof2_table']] = gof2_table 
			allData[['gof_sfs']] = gof_sfs 
			allData[['gof2_sfs']] = gof2_sfs 
			allData[['meta']] = meta
			
			return(allData)
		}
	})

	## observed sum stats, demographic inferences
	output$visualization_data <- renderUI({
		if(is.null(input$results) == FALSE){
			tabsetPanel(
				type = "tabs",
				tabPanel("Observed summary statistics", uiOutput("user_dataset_tabset")),
				tabPanel("Demographic inferences", uiOutput("user_inferences")),
				tabPanel("Definitions of statistics and parameters", uiOutput("definitions"))
			)
		}else{
			return()
		}
	})

	output$user_dataset_tabset <- renderUI({
		if(is.null(input$results) == FALSE){
			if(allData()[['users_infos']][1,2]==2){
				# if number of species == 2
				tabsetPanel(id = "observed_dataset",
				type = "tabs",
				tabPanel("Summarized jSFS",
					fluidRow(
						column( width = 12, style='padding:20px;',
							prettyCheckbox(inputId = "show_points_stats_sites", shape = "round", value = FALSE, label = strong("Show individual loci"), icon = icon("check"), animation = "tada", status = "success", bigger = TRUE)
							),
						
						column( width = 12, style="margin-top:-0.5em",
							plotlyOutput("plot_obs_stats_sites")
							)
						)
				),
				
				tabPanel("Polymorphism",
					fluidRow(
						column( width = 12, style='padding:20px;',
							prettyCheckbox(inputId = "show_points_diversity", shape = "round", value = FALSE, label = strong("Show individual loci"), icon = icon("check"), animation = "tada", status = "success", bigger = TRUE)
							),
						
						column( width = 12, style="margin-top:-0.5em",
							plotlyOutput("plot_obs_stats_diversity")
							)
						)
				),

				tabPanel("Tajima's D",
					fluidRow(
						column( width = 12, style='padding:20px;',
							prettyCheckbox(inputId = "show_points_Tajima", shape = "round", value = FALSE, label = strong("Show individual loci"), icon = icon("check"), animation = "tada", status = "success", bigger = TRUE)
							),
						
						column( width = 12, style="margin-top:-0.5em",
							plotlyOutput("plot_obs_stats_tajima")
							)
						)
				),
				
				tabPanel("Differentiation and divergence",
					fluidRow(
						column( width = 12, style='padding:20px;',
							prettyCheckbox(inputId = "show_points_divergence", shape = "round", value = FALSE, label = strong("Show individual loci"), icon = icon("check"), animation = "tada", status = "success", bigger = TRUE)
							),
						
						column( width = 12, style="margin-top:-0.5em",
							plotlyOutput("plot_obs_stats_divergence")
							)
						)
					)
				)
			}else{
				if(allData()[['users_infos']][1,2]==1){
					# if number of species == 1
					fluidRow(
						column( width = 12, style='padding:20px;',
							prettyCheckbox(inputId = "show_points_stats_sites_1pop", shape = "round", value = FALSE, label = strong("Show individual loci"), icon = icon("check"), animation = "tada", status = "success", bigger = TRUE)
							),
						
						column( width = 12, style="margin-top:-0.5em",
							plotlyOutput("plot_obs_stats_1pop")
							)
						)
				}else{
					return()
				}
			}
		}else{
			return()
		}
	})
	
	output$user_inferences <- renderUI({
		if(is.null(input$results) == FALSE){
			if(allData()[['users_infos']][1,2]==2){
				# if number of species == 2
				tabsetPanel(id = "inferences",
				type = "tabs",
				tabPanel("Multilocus model comparison", uiOutput("display_modComp")),
				tabPanel("Locus specific model comparison", numericInput("threshold_locus_specific_model_comp", label = h3("Posterior probability threshold value below which an inference is considered ambiguous"), width = (0.25*as.numeric(input$dimension[1])), value = 0.9, min = 0, max = 1, step = 0.005), hr(), plotlyOutput("locus_specific_model_comparison", height = 'auto', width = 'auto')),
				tabPanel("Estimated parameters", uiOutput("parameters_estimates")),
				tabPanel("Goodness-of-fit test", uiOutput("gof"))
				)
			}else{
				if(allData()[['users_infos']][1,2]==1){
					# if number of species == 1
					tabsetPanel(id = "inferences",
					type = "tabs",
					tabPanel("Multilocus model comparison", uiOutput("display_modComp")),
					tabPanel("Estimated parameters", uiOutput("parameters_estimates")),
					tabPanel("Goodness-of-fit test", uiOutput("gof"))
					)
				}
			}
		}else{
			return()
		}
	})

	output$definitions <- renderUI({
		if(is.null(input$results) == FALSE){
				# if number of species == 2
				tabsetPanel(id = "definitions",
					type = "tabs",
					tabPanel('Summary statistics', uiOutput('definitions_statistics')),
					tabPanel('Model parameters', uiOutput('definitions_parameters'))
				)
		}
	})
	
	output$definitions_parameters <- renderUI({
		if(is.null(input$results)){
			return(NULL)
		}else{
			fluidPage(
				HTML('<h2><u>Population sizes</u></h2>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Na</b>: effective size of the ancestral population [# of diploid individuals]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>N1</b>; <b>N2</b>: effective size of population 1 (resp. 2) [# of diploid individuals]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<u>If <i>Ne</i> is genomically heterogeneous</u></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>shape_N_a</b>; <b>shape_N_b</b>: shape parameter &alpha; (resp. &beta;) of the Beta(&alpha;, &beta;) distribution for <i>Ne</i> (shared by all populations)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<u>If there is a demographic change</u></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Tdem1</b>; <b>Tdem2</b>: time of the demographic change in population 1 (resp. 2) [# of generations]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>founders1</b>; <b>founders2</b>: number of founder individuals in population 1 (resp. 2) at the time of the demographic change <i>T<sub>dem1</sub></i> (and <i>T<sub>dem2</sub></i>)</h3>'),

				HTML('<h2><u>Times for two populations or more</u></h2>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Tsplit</b>: time of split at which the ancestral population subdivides in two populations [# of generations]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<u>If there is both gene flow and isolation</u></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Tsc</b>: time of secondary contact at which the two populations start exchanging genes [# of generations]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>Tam</b>: time of ancient migration at which the two populations stop exchanging genes [# of generations]</h3>'),

				HTML('<h2><u>Migration and barriers</u></h2>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>M12</b>; <b>M21</b>: introgression rate from population 2 to 1 (resp. from 1 to 2) [# of migrants per generation]</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<u>If <i>N.m</i> is genomically heterogeneous (bimodal model)</u></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>nBarriersM12</b>; <b>nBarriersM21</b>: number of loci inferred as interspecies barriers for introgression from population 2 to 1 (resp. from 1 to 2)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<u>If <i>N.m</i> is genomically heterogeneous (beta model)</u></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>shape_M12_a</b>; <b>shape_M12_b</b>: shape parameter &alpha; (resp. &beta;) of the Beta(&alpha;, &beta;) distribution for <i>N.m</i> (from population 2 to 1)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>shape_M21_a</b>; <b>shape_M21_b</b>: shape parameter &alpha; (resp. &beta;) of the Beta(&alpha; &beta;) distribution for <i>N.m</i> (from population 1 to 2)</h3>'),
		)}
	})

	output$definitions_statistics <- renderUI({
		if(is.null(input$results)){
			return(NULL)
		}else{
			fluidPage(
				HTML('<h2><u>Summary statistics</u></h2>'),
				HTML('<h3><b>dataset</b>: name of the target locus</h3>'),
				HTML('<h3><b><i>Summarized jSFS</i></b></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>sf_avg</b>: fraction of sites with a fixed difference between the populations/species</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>sxA_avg</b>; <b>sxB_avg</b>: fraction of sites with a polymorphism specific to each population/species</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>ss_avg</b>: fraction of sites with a polymorphism shared between the population/species</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>successive_ss_avg</b>: maximal number of successive shared sites in the target locus</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>ss_sf</b>: if the target locus has at least one shared site (ss) and one fixed difference (sf), the value is set to 1; and 0 otherwise</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>ss_noSf</b>: if the target locus has at least one shared site (ss) but no fixed difference (sf), the value is set to 1; and 0 otherwise</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>noSs_sf</b>: if the target locus has no shared site (ss) but at least one fixed difference (sf), the value is set to 1; and 0 otherwise</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>noSs_noSf</b>: if the target locus has no shared site (ss) and no fixed difference (sf), the value is set to 1; and 0 otherwise</h3>'),
				HTML('<h3><b><i>Polymorphism</i></b></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>piA_avg</b>; <b>piB_avg</b>: pairwise nucleotide diversity (&pi;) for each population/species (Tajima, 1983)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>thetaA_avg</b>; <b>thetaB_avg</b>: Watterson’s theta for each population/species (Watterson, 1975)</h3>'),
				HTML('<h3><b><i>Tajima&#39;s D</i></b></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>DtajA_avg</b>; <b>DtajB_avg</b>: Tajima&#39;s <i>D</i> for each population/species (Tajima, 1989)</h3>'),
				HTML('<h3><b><i>Differentiation and divergence</i></b></h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>divAB_avg</b>: raw divergence (<i>D<sub>xy</sub></i>) between the population/species (Nei, 1987)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>netdivAB_avg</b>: net divergence (<i>D<sub>a</sub></i>) between the populations/species, measured by <i>D<sub>xy</sub></i> - (&pi;<sub>A</sub> + &pi;<sub>B</sub>)/2 (Nei & Li, 1979)</h3>'),
				HTML('<h3>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>FST_avg</b>: <i>F<sub>ST</sub></i> measured by 1-&pi;<sub>S</sub>/&pi;<sub>T</sub>; where &pi;<sub>S</sub> is the average nucleotide diversity in each population/species and &pi;<sub>T</sub> is the total nucleotide diversity over the populations/species (Wright, 1943)</h3>'),
				hr(),
				HTML('<h2><u>References</u></h2>'),
				HTML('<h3>Nei, M. (1987). Molecular Evolutionary Genetics. Columbia University Press, New York.'),
				HTML('<h3>Nei, M. & Li, W‐H. (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. PNAS, 76: 5269–5273.'),
				HTML('<h3>Tajima, F. (1989). The effect of change in population size on DNA polymorphism. Genetics, 123(3): 597-601.'),
				HTML('<h3>Tajima, F. (1983). Evolutionary relationship of DNA sequences in finite populations. Genetics, 105(2): 437-460.'),
				HTML('<h3>Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. Theor. Popul. Biol., 7(2): 256-276.'),
				HTML('<h3>Wright, S. (1943). Isolation by distance. Genetics, 28: 114–138.')

			)
		}
	})
	
	# Display the model comparisons	
	output$display_modComp <- renderUI({
		if(is.null(allData()[['hierarchical']])){
			return(NULL)
		}else if(allData()[['users_infos']][1,2] == 2){
			# if number of species == 2
			 if(allData()[['hierarchical']][2,1] == 'isolation'){
				fluidPage(
					hr(),
					infoBox("Migration versus isolation", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[1], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[1], sep=''), icon = icon("check"), color='navy'),
					infoBox("AM versus SI", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[2], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[2], sep=''), icon = icon("check"), color='navy'),
					infoBox("N-homo versus N-hetero", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[3], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[3], sep=''), icon = icon("check"), color='navy')
				)
			}else if(allData()[['hierarchical']][2,1] == 'migration'){
				fluidPage(
					hr(),
					infoBox("Migration versus isolation", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[1], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[1], sep=''), icon = icon("check"), color='navy'),
					infoBox("IM versus SC", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[2], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[2], sep=''), icon = icon("check"), color='navy'),
					infoBox("N-homo versus N-hetero", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[3], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[3], sep=''), icon = icon("check"), color='navy'),
					infoBox("M-homo versus M-hetero", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[4], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[4], sep=''), icon = icon("check"), color='navy')
				)
			}
		}
		else if(allData()[['users_infos']][1,2] == 1){
			# if number of species == 1
			fluidPage(
				hr(),
				infoBox("Expansion versus Constant versus Contraction", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[1], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[1], sep=''), icon = icon("check"), color='navy'),
				infoBox("N-homo versus N-hetero", paste('best model = ', as.matrix(allData()[['hierarchical']][2,])[2], sep=''), paste('post. proba = ', round(as.numeric(as.matrix(allData()[['hierarchical']][3,])), 5)[2], sep=''), icon = icon("check"), color='navy')
			)
		}
	})


	output$prior_mutation <- renderUI({
		fluidPage(
			boxPlus(title = h2("Mutation and recombination"), height = NULL, width = NULL, closable = FALSE, status = "primary", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				fluidRow(
					column(width=5, numericInput("mu", label = h5('Mutation rate'), value = 0.000000003, min = 0, max = 0.00001)),
					column(width=5, numericInput("rho_over_theta", label = h5('Ratio r/µ'), value = 0.1, min = 0, max = 5))
				)
			),
			
			boxPlus(
				title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "teal",
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "INFORMATION", label_status = "primary",
				
				h3("The mutation rate ", strong("(µ)"), "is", strong("the probability per generation and per nucleotide"), "that an allele will not be properly replicated."),
				h3("If an external group ", strong("is not specified"), "then all genes/contigs/locus share the same µ."),
				HTML("<h3>If an external group is specified then the local µ, for a locus <i>i</i> is corrected by µ * <i>div<sub>i</sub></i> / <i>div<sub>avg</sub></i> where <i>div<sub>i</sub></i> is the local divergence between the ingroup and the outgroup at that locus, and <i>div<sub>avg</sub></i> is the divergence averaged over all loci</h3>"),
				br(),
				h3("The r/µ ratio is the ratio of recombination (/bp /generation) over mutation (/bp /generation).")
			)
		)
	})


	output$prior_Ne <- renderUI({
		fluidPage(
			boxPlus(title = h2("Population size"), height = 300, width = NULL, closable = FALSE, status = "danger", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
				fluidRow(
					column(width=5, numericInput("N_min", label = h5('min'), value = 100, min = 0, max = 19999999)),
					column(width=5, numericInput("N_max", label = h5('max'), value = 1000000, min = 1, max = 20000000))
				)
			),
			
			boxPlus(
				title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "danger",
				boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
				enable_label = TRUE, label_text = "INFORMATION", label_status = "danger",
				
				h3("The effective population size ", em(strong("Ne")), "is the number of diploid individuals within current and ancestral species/populations."),
				hr(),
				h3("In the", strong("ABC"), "simulations,", em(strong("Ne")), "will be drawn from the setted prior distribution independently for all current and ancestral species/populations")
			)

		)
	})


	output$prior_times <- renderUI({
		if(input$nspecies==2){
			fluidPage(
				boxPlus(title = h2("Time of split"), height = NULL, width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
					fluidRow(
						column(width=5, numericInput("Tsplit_min", label = h5('min'), value = 100, min = 0, max = 49999999)),
						column(width=5, numericInput("Tsplit_max", label = h5('max'), value = 1000000, min = 1, max = 50000000))
					)
				),

				boxPlus(
					title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "warning",
					boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
					enable_label = TRUE, label_text = "INFORMATION", label_status = "warning",
					
					HTML("<h3>The speciation time <b>T<sub>split</sub></b> is expressed <b>in number of generations.</b></h3>"),
					h3("For annual organisms: one generation = one year."),
					h3("For perennial organisms: one generation = average age for an individual to transmit a descendant (which is different from the age of sexual maturity)."),
					br(),
					HTML("<h3>The prior distribution is uniform between <b>T<sub>split_min</sub></b> and <b>T<sub>split_max</sub></b>.</h3>"),
					HTML("<h3>For each simulation in the <b>SC</b> and <b>AM</b> models, the time of secondary contact between lines (T<sub>SC</sub>) and and arrest of ancient migration (T<sub>AM</sub>) are drawn uniformly between <b>T<sub>split_min</sub></b> and the sampled T<sub>split</sub></h3>")
				)
			)
		}else{
			fluidPage(
				boxPlus(title = h2("Time of demographic change"), height = NULL, width = NULL, closable = FALSE, status = "warning", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
					fluidRow(
						column(width=5, numericInput("Tsplit_min", label = h5('min'), value = 100)),
						column(width=5, numericInput("Tsplit_max", label = h5('max'), value = 1000000))
					)
				),

				boxPlus(
					title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = "warning",
					boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
					enable_label = TRUE, label_text = "INFORMATION", label_status = "warning",
					
					h3("The time of demographic change", strong(em("Tdem")), "is expressed", strong("in number of generations.")),
					h3("For annual organisms: one generation = one year."),
					h3("For perennial organisms: one generation = average age for an individual to transmit a descendant (which is different from the age of sexual maturity)."),
					h3("It represents the number of generations since population expansion or contraction.")
				)
			)
		}
	})

	output$prior_migration <- renderUI({
		if(input$nspecies==2){
			fluidPage(
				boxPlus(title = h2("Migration rates"), height = 300, width = NULL, closable = FALSE, status = "success", solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
					fluidRow(
						column(width=5, numericInput("M_min", label = h5('min'), value = 0.4)),
						column(width=5, numericInput("M_max", label = h5('max'), value = 20))
					),
					fluidRow(
						column(width=5,	selectInput("modeBarrier", label = h4("Model for barriers"), choices = list("bimodal" = 'bimodal', "beta" = 'beta'), selected = 'bimodal'))
					)
				),
			
				boxPlus(
					title = "", width = NULL, icon = NULL, solidHeader = TRUE, gradientColor = NULL,
					boxToolSize = "lg", footer_padding = TRUE, collapsible = TRUE, collapsed = TRUE, closable = FALSE,
					enable_label = TRUE, label_text = "INFORMATION", label_status = "success",
					HTML("<h3>Migration rates are expressed in <b>4.<i>Ne.m</i></b> where <b><i>m</i></b> is the fraction of each subpopulation made up of new migrants each generation.</h3>")
				)
			)
		}
	})


	output$output_posterior_2pops <- renderUI({
		fileName = input$results
		if(is.null(fileName)) {
			return()
		}else{
			fluidPage(
				fluidRow( width = 12,
					selectInput('param_name', 'Parameter: ', allData()[['list_parameters']])
				),
				
				fluidRow( width = 12,
					plotlyOutput( outputId = 'posterior_parameters_2pops')
				)
			)
		}
	})

	output$output_posterior_1pop <- renderUI({
		fileName = input$results
		if(is.null(fileName)) {
			return()
		}else{
			fluidPage(
				fluidRow( width = 12,
					selectInput('param_name', 'Parameter: ', allData()[['list_parameters']])
				),
				
				fluidRow( width = 12,
					plotlyOutput( outputId = 'posterior_parameters_1pop')
				)
			)
		}
	})

	
	output$distribution_gof <- renderUI({
		fileName = input$results
		if(is.null(fileName)) {
			return()
		}else{
			fluidPage(
				fluidRow( width = 12,
					selectInput('stat_name', 'Summary statistic: ', list_statistics())
				),
				
				fluidRow( width = 12,
					plotlyOutput( outputId = 'density_statistic')
				)
			)
		}
	})


	## Get the names of the parameters
	list_statistics = reactive({
		fileName = input$results
		if(is.null(fileName)){
			return (NULL)
		}else{
			rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
			x = allData()[['distribution_PCA']] 
			
			toRemove = c(1)
			for(i in 1:(ncol(x)-1)){ # -1 to keep the column prior/posterior/optimized/observed
				if(sd(x[,i]) < 0.00001){
					toRemove = c(toRemove, i)
				}
			}
			toRemove = c(toRemove, grep('max', colnames(x)))
			toRemove = c(toRemove, grep('min', colnames(x)))
			toRemove = c(toRemove, grep('origin', colnames(x)))
			toRemove = unique(toRemove)
			x = x[, -toRemove]

			# remove the temporary unarchived results
			return(colnames(x))
		}
	})


	##  Plot the posterior
	output$posterior_parameters_2pops <- renderPlotly({
		fileName = input$results
		if (is.null(fileName)){
			return(NULL)
		}else{
			param_name = input$param_name
			
			rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]

			theme_set(theme_classic())
			figure = list()

#			bpfile = read.table(paste(rootName, '/bpfile', sep=''), h=F, skip=1); nLoci = ncol(bpfile)
			yaml = allData()[['yaml']]
			Nref = allData()[['Nref']]
			priorfile = allData()[['priorfile']]
			res1 = allData()[['posterior']]
			res2 = allData()[['optimized_posterior']]
			
			# prior
			if( param_name == 'N1' ){
				scale = Nref
				prior = priorfile$N1*Nref
			}
			if( param_name == 'N2' ){
				scale = Nref
				prior = priorfile$N2*Nref
			}
			
			if( param_name == 'Na' ){
				scale = Nref
				prior = priorfile$Na*Nref
			}
		
			if( param_name == 'founders1' ){
				scale = 1
				prior = priorfile$founders1 
			}

			if( param_name == 'founders2' ){
				scale = 1
				prior = priorfile$founders2 
			}

			if( param_name == 'M12'){
				scale = 1
				prior = priorfile$M12
			}
			
			if( param_name == 'M21'){
				scale = 1
				prior = priorfile$M21
			}

			if( param_name == 'Tsplit' || param_name == 'Tam' || param_name == 'Tsc' || param_name == 'Tdem1' || param_name == 'Tdem2' ){
				scale = 4*Nref
				
				if( param_name == 'Tsplit' ){
					prior = priorfile$Tsplit*scale
				}

				if( param_name == 'Tam' ){
					prior = priorfile$Tam*scale
				}
				
				if( param_name == 'Tsc' ){
					prior = priorfile$Tsc*scale
				}
				
				if( param_name == 'Tdem1'){
					prior = priorfile$Tdem1*scale
				}
				
				if( param_name == 'Tdem2'){
					prior = priorfile$Tdem2*scale
				}
			}

			if( param_name == "shape_N_a" || param_name == "shape_N_b" || param_name == "shape_M12_a" || param_name == "shape_M12_b" || param_name == "shape_M21_a" || param_name == "shape_M21_b" ){
				scale = 1
				if( param_name == 'shape_N_a'){
					prior = priorfile$shape_N_a
				}
				
				if( param_name == 'shape_N_b'){
					prior = priorfile$shape_N_b
				}
				
				if( param_name == 'shape_M12_a'){
					prior = priorfile$shape_M12_a
				}
				
				if( param_name == 'shape_M12_b'){
					prior = priorfile$shape_M12_b
				}
				
				if( param_name == 'shape_M21_a'){
					prior = priorfile$shape_M12_a
				}
				
				if( param_name == 'shape_M21_b'){
					prior = priorfile$shape_M12_b
				}
			}

			if( param_name == "nBarriersM12" ){
				scale = 1
				prior = priorfile$nBarriersM12
			}
			
			if( param_name == "nBarriersM21" ){
				scale = 1
				prior = priorfile$nBarriersM21
			}
			
			#	prior = x[,i] * scale
			#	prior = data.frame(x = prior, label=rep('prior', length(prior)))
			prior = data.frame(x = prior, distribution=rep("Prior", length(prior)))
			posterior1 = res1[,which(colnames(res1)==param_name)] * scale
			posterior1 = data.frame(x = posterior1, distribution=rep('Posterior', length(posterior1)))
			
			posterior2 = res2[,which(colnames(res1)==param_name)] * scale
			posterior2 = data.frame(x = posterior2, distribution=rep('Optimized posterior', length(posterior2)))
#			
			df=rbind(prior, posterior1, posterior2)


			p <- ggplot(df, aes(x, fill = distribution)) + geom_density(alpha = 0.7, size = 0.25) + scale_fill_manual(values=c("darkgray", viridis_pal(option="D")(5)[c(1,4)])) + theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), legend.text = element_text(size = 15)) + scale_x_continuous(name = param_name)
			p <- ggplotly(p, width = 0.55*as.numeric(input$dimension[1]), height = 0.55*as.numeric(input$dimension[2]))

			figure_estimations = ggplotly(p)
			return(figure_estimations)
		}
	})

	# posterior 1 pop
	##  Plot the posterior
	output$posterior_parameters_1pop <- renderPlotly({
		fileName = input$results
		if (is.null(fileName)){
			return(NULL)
		}else{
			param_name = input$param_name
		
			theme_set(theme_classic())
			priorfile = allData()[['priorfile']]
			res1 = allData()[['posterior']]
			res2 = allData()[['optimized_posterior']]

			prior = priorfile[,which(colnames(priorfile)==param_name)]
			prior = data.frame(x = prior, distribution=rep('Prior', length(prior)))
			
			posterior1 = res1[,which(colnames(res1)==param_name)]
			posterior1 = data.frame(x = posterior1, distribution=rep('Posterior', length(posterior1)))
			
			posterior2 = res2[,which(colnames(res2)==param_name)]
			posterior2 = data.frame(x = posterior2, distribution=rep('Optimized posterior', length(posterior2)))
			
			df=rbind(prior, posterior1, posterior2)

			p <- ggplot(df, aes(x, fill = distribution)) + geom_density(alpha = 0.7, size = 0.25) + scale_fill_manual(values=c("white", viridis_pal(option="D")(5)[c(1,4)])) + theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), legend.text = element_text(size = 15)) + scale_x_continuous(name = param_name)
			p <- ggplotly(p, width = 0.55*as.numeric(input$dimension[1]), height = 0.55*as.numeric(input$dimension[2]))

			figure_estimations = ggplotly(p)
			return(figure_estimations)
		}
	})


	# Display the table with estimated parameters	
	output$table_parameters_2pops <- renderPlotly({
		fileName = input$results
		if (is.null(fileName)){
			return(NULL)
		}else{
			rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
		
			# table
			theme_set(theme_classic())
			figure = list()

			# read information
			Nref = allData()[['Nref']]
			if(input$modeInference == 'neural network'){
				res1 = allData()[['posterior']]
				res2 = allData()[['optimized_posterior']]
				title_plot = 'ABC - neural network'
				# table
				param_names = NULL
				post1_Q1 = NULL
				post1_median = NULL
				post1_Q2 = NULL
				post2_Q1 = NULL
				post2_median = NULL
				post2_Q2 = NULL

				nparams = ncol(res1)
				for(i in 1:nparams){
					param_name = colnames(res1)[i]
					
					scale = 1
					if( param_name == 'N1' || param_name == 'N2' || param_name == 'Na' ){
						scale = Nref
					}
					
					if( param_name == 'Tsplit' || param_name == 'Tam' || param_name == 'Tsc' || param_name == 'Tdem1' || param_name == 'Tdem2' ){
						scale = 4*Nref
					}
					
					if( param_name == 'M12' || param_name == 'M21'){
						scale = 1/4
					}
						
					posterior1 = res1[,i] * scale
					posterior1 = data.frame(x = posterior1, label=rep('Posterior', length(posterior1)))

					posterior2 = res2[,i] * scale
					posterior2 = data.frame(x = posterior2, label=rep('First optimization', length(posterior2)))
					
					# table
					param_names = c(param_names, param_name)
					if(param_name%in%c('N1', 'N2', 'Na', 'Tsplit', 'Tsc', 'Tam', 'Tmin', 'Tdem1', 'Tdem2', 'nBarriersM12', 'nBarriersM21')){
						post1_Q1_tmp = formatC(round(quantile(posterior1$x, 0.025), 0), format='d', big.mark=' ')
						post1_median_tmp = formatC(round(quantile(posterior1$x, 0.5), 0), format='d', big.mark=' ')
						post1_Q2_tmp = formatC(round(quantile(posterior1$x, 0.975), 0), format='d', big.mark=' ')

						post2_Q1_tmp = formatC(round(quantile(posterior2$x, 0.025), 0), format='d', big.mark=' ')
						post2_median_tmp = formatC(round(quantile(posterior2$x, 0.5), 0), format='d', big.mark=' ')
						post2_Q2_tmp = formatC(round(quantile(posterior2$x, 0.975), 0), format='d', big.mark=' ')
					}else{
						post1_Q1_tmp = formatC(round(quantile(posterior1$x, 0.025), 5), format="f", big.mark=" ", digits=5)
						post1_median_tmp = formatC(round(quantile(posterior1$x, 0.5), 5), format="f", big.mark=" ", digits=5)
						post1_Q2_tmp = formatC(round(quantile(posterior1$x, 0.975), 5), format="f", big.mark=" ", digits=5)
						
						post2_Q1_tmp = formatC(round(quantile(posterior2$x, 0.025), 5), format="f", big.mark=" ", digits=5)
						post2_median_tmp = formatC(round(quantile(posterior2$x, 0.5), 5), format="f", big.mark=" ", digits=5)
						post2_Q2_tmp = formatC(round(quantile(posterior2$x, 0.975), 5), format="f", big.mark=" ", digits=5)
						
					}
						
					post1_Q1 = c(post1_Q1, post1_Q1_tmp)
					post1_median = c(post1_median, post1_median_tmp)
					post1_Q2 = c(post1_Q2, post1_Q2_tmp)
					
					post2_Q1 = c(post2_Q1, post2_Q1_tmp)
					post2_median = c(post2_median, post2_median_tmp)
					post2_Q2 = c(post2_Q2, post2_Q2_tmp)
				}

				# print table
				col_tmp = viridis_pal(option="D", alpha=1)(5)
				col_post1_header = col_tmp[1]
				col_post2_header = col_tmp[4]
				col_tmp = viridis_pal(option="D", alpha=0.4)(5)
				col_post1 = col_tmp[1]
				col_post2 = col_tmp[4]
				green = "#C7F464"
				dark_grey = "#1e2b37"
				light_grey = "#556270"

				table_estimations = plot_ly( type = 'table',
					header = list(
						values = c("<b>Parameter</b>", "<b>HPD 0.025</b>", "<b>HPD median (posterior)</b>", "<b>HPD 0.975</b>", "<b>HPD 0.025</b>", "<b>HPD median (optimized posterior)</b>", "<b>HPD 0.975</b>"),
						line = list(color = dark_grey),
						fill = list(color = c(dark_grey, col_post1_header, col_post1_header, col_post1_header, col_post2_header, col_post2_header, col_post2_header)),
						align = c('left','center'),
						font = list(color = c(green, rep("white", 6), size = 30))
					),
					cells = list(
						values = rbind(
							paste('<b>', param_names, '</b>', sep=''),
							post1_Q1,
							paste('<b>', post1_median, '</b>', sep=''),
							post1_Q2,
							
							post2_Q1,
							paste('<b>', post2_median, '</b>', sep=''),
							post2_Q2
						),
						line = list(color = dark_grey),
						fill = list(color = c(light_grey, col_post1, col_post1, col_post1, col_post2, col_post2, col_post2)),
						align = c('left', 'center'),
						font = list(color = c(green, dark_grey,  size = 30))
					), 
					width = 0.75*as.numeric(input$dimension[1]), height = 0.75*as.numeric(input$dimension[2])
				) %>% layout(title = title_plot) %>%
					layout(plot_bgcolor='rgb(1,1,1,0)') %>% 
					layout(paper_bgcolor='rgb(1,1,1,0)') 
				return(table_estimations)
			}else{
				res1 = allData()[['posterior_RF']]
				res2 = allData()[['optimized_posterior_RF']]
				title_plot = 'ABC - random forest'
				# table
				param_names = NULL
				post1_Q1 = NULL
				post1_median = NULL
				post1_Q2 = NULL
				post2_Q1 = NULL
				post2_median = NULL
				post2_Q2 = NULL

				nparams = nrow(res1)
				for(i in 1:nparams){
					param_name = res1[i,1]
					
					scale = 1
					if( param_name == 'N1' || param_name == 'N2' || param_name == 'Na' ){
						scale = Nref
					}
					
					if( param_name == 'Tsplit' || param_name == 'Tam' || param_name == 'Tsc' || param_name == 'Tdem1' || param_name == 'Tdem2' ){
						scale = 4*Nref
					}
					
					if( param_name == 'M12' || param_name == 'M21'){
						scale = 1/4
					}
						
					posterior1 = as.numeric(res1[i,-1]) * scale
					posterior2 = as.numeric(res2[i,-1]) * scale
					
					# table
					param_names = c(param_names, as.character(param_name))
					if(param_name%in%c('N1', 'N2', 'Na', 'Tsplit', 'Tsc', 'Tam', 'Tmin', 'Tdem1', 'Tdem2', 'nBarriersM12', 'nBarriersM21')){
						post1_Q1_tmp = formatC(round(posterior1[1], 0), format='d', big.mark=' ')
						post1_median_tmp = formatC(round(posterior1[2], 0), format='d', big.mark=' ')
						post1_Q2_tmp = formatC(round(posterior1[3], 0), format='d', big.mark=' ')

						post2_Q1_tmp = formatC(round(posterior2[1], 0), format='d', big.mark=' ')
						post2_median_tmp = formatC(round(posterior2[2], 0), format='d', big.mark=' ')
						post2_Q2_tmp = formatC(round(posterior2[3], 0), format='d', big.mark=' ')
					}else{
						post1_Q1_tmp = formatC(round(posterior1[1], 5), format="f", big.mark=" ", digits=5)
						post1_median_tmp = formatC(round(posterior1[2], 5), format="f", big.mark=" ", digits=5)
						post1_Q2_tmp = formatC(round(posterior1[3], 5), format="f", big.mark=" ", digits=5)
						
						post2_Q1_tmp = formatC(round(posterior2[1], 5), format="f", big.mark=" ", digits=5)
						post2_median_tmp = formatC(round(posterior2[2], 5), format="f", big.mark=" ", digits=5)
						post2_Q2_tmp = formatC(round(posterior2[3], 5), format="f", big.mark=" ", digits=5)
					}
						
					post1_Q1 = c(post1_Q1, post1_Q1_tmp)
					post1_median = c(post1_median, post1_median_tmp)
					post1_Q2 = c(post1_Q2, post1_Q2_tmp)
					
					post2_Q1 = c(post2_Q1, post2_Q1_tmp)
					post2_median = c(post2_median, post2_median_tmp)
					post2_Q2 = c(post2_Q2, post2_Q2_tmp)
				}

				# print table
				col_tmp = viridis_pal(option="D", alpha=1)(5)
				col_post1_header = col_tmp[1]
				col_post2_header = col_tmp[4]
				col_tmp = viridis_pal(option="D", alpha=0.4)(5)
				col_post1 = col_tmp[1]
				col_post2 = col_tmp[4]
				green = "#C7F464"
				dark_grey = "#1e2b37"
				light_grey = "#556270"

				table_estimations = plot_ly( type = 'table',
					header = list(
						values = c("<b>Parameter</b>", "<b>HPD 0.025</b>", "<b>HPD median (posterior)</b>", "<b>HPD 0.975</b>", "<b>HPD 0.025</b>", "<b>HPD median (optimized posterior)</b>", "<b>HPD 0.975</b>"),
						line = list(color = dark_grey),
						fill = list(color = c(dark_grey, col_post1_header, col_post1_header, col_post1_header, col_post2_header, col_post2_header, col_post2_header)),
						align = c('left','center'),
						font = list(color = c(green, rep("white", 6), size = 30))
					),
					cells = list(
						values = rbind(
							paste('<b>', param_names, '</b>', sep=''),
							post1_Q1,
							paste('<b>', post1_median, '</b>', sep=''),
							post1_Q2,
							
							post2_Q1,
							paste('<b>', post2_median, '</b>', sep=''),
							post2_Q2
						),
						line = list(color = dark_grey),
						fill = list(color = c(light_grey, col_post1, col_post1, col_post1, col_post2, col_post2, col_post2)),
						align = c('left', 'center'),
						font = list(color = c(green, dark_grey,  size = 30))
					), 
					width = 0.75*as.numeric(input$dimension[1]), height = 0.75*as.numeric(input$dimension[2])
				) %>% layout(title = "ABC - Random Forest") %>% 
					layout(plot_bgcolor='rgb(1,1,1,0)') %>% 
					layout(paper_bgcolor='rgb(1,1,1,0)') 
				return(table_estimations)
			}
			
		}
	})
	
	output$table_parameters_1pop <- renderPlotly({
		fileName = input$results
		if (is.null(fileName)){
			return(NULL)
		}else{
			# table
			theme_set(theme_classic())
			figure = list()

			# read information
			modeInference = input$modeInference # RF or nnet
			if(modeInference == 'neural network'){
				retained_values = seq(1, nrow(allData()[['table_optimized_posterior']]), 2)
				title_plot = 'ABC - neural network'
			}else{
				retained_values = seq(2, nrow(allData()[['table_optimized_posterior']]), 2)
				title_plot = 'ABC - random forest'
			}
			
			res1 = allData()[['table_posterior']][retained_values,]
			res4 = allData()[['table_optimized_posterior']][retained_values,]
			
			# table
			param_names = NULL
			post1_Q1 = NULL
			post1_median = NULL
			post1_Q2 = NULL
			
			post4_Q1 = NULL
			post4_median = NULL
			post4_Q2 = NULL

			nparams = nrow(res1)
			
			for(i in 1:nparams){
				param_name = as.character(res1[i,1])
				
				# table
				param_names = c(param_names, param_name)
				if(param_name%in%c('N', 'Npast', 'Tdem')){
					post1_Q1_tmp = formatC(round(res1[i,2], 0), format='d', big.mark=' ')
					post1_median_tmp = formatC(round(res1[i,3], 0), format='d', big.mark=' ')
					post1_Q2_tmp = formatC(round(res1[i,4], 0), format='d', big.mark=' ')
					
					post4_Q1_tmp = formatC(round(res4[i,2], 0), format='d', big.mark=' ')
					post4_median_tmp = formatC(round(res4[i,3], 0), format='d', big.mark=' ')
					post4_Q2_tmp = formatC(round(res4[i,4], 0), format='d', big.mark=' ')
				}else{
					post1_Q1_tmp = formatC(round(res1[i,2], 5), format='f', big.mark=' ')
					post1_median_tmp = formatC(round(res1[i,3], 5), format='f', big.mark=' ')
					post1_Q2_tmp = formatC(round(res1[i,4], 5), format='f', big.mark=' ')
					
					post4_Q1_tmp = formatC(round(res4[i,2], 5), format='f', big.mark=' ')
					post4_median_tmp = formatC(round(res4[i,3], 5), format='f', big.mark=' ')
					post4_Q2_tmp = formatC(round(res4[i,4], 5), format='f', big.mark=' ')
				}
					
				post1_Q1 = c(post1_Q1, post1_Q1_tmp)
				post1_median = c(post1_median, post1_median_tmp)
				post1_Q2 = c(post1_Q2, post1_Q2_tmp)
				
				post4_Q1 = c(post4_Q1, post4_Q1_tmp)
				post4_median = c(post4_median, post4_median_tmp)
				post4_Q2 = c(post4_Q2, post4_Q2_tmp)
			}

			# print table
			col_tmp = viridis_pal(option="D", alpha=1)(5)
			col_post1_header = col_tmp[1]
			col_post4_header = col_tmp[4]
			col_tmp = viridis_pal(option="D", alpha=0.4)(5)
			col_post1 = col_tmp[1]
			col_post4 = col_tmp[4]
			green = "#C7F464"
			dark_grey = "#1e2b37"
			light_grey = "#556270"

			table_estimations = plot_ly( type = 'table',
				header = list(
					values = c("<b>Parameter</b>", "<b>HPD 0.025</b>", "<b>HPD median (posterior)</b>", "<b>HPD 0.975</b>", "<b>HPD 0.025</b>", "<b>HPD median (optimized posterior)</b>", "<b>HPD 0.975</b>"),
					line = list(color = dark_grey),
					fill = list(color = c(dark_grey, col_post1_header, col_post1_header, col_post1_header, col_post4_header, col_post4_header, col_post4_header)),
					align = c('left','center'),
					font = list(color = c(green, rep("white", 6), size = 30))
				),
				cells = list(
					values = rbind(
						paste('<b>', param_names, '</b>', sep=''),
						post1_Q1,
						paste('<b>', post1_median, '</b>', sep=''),
						post1_Q2,
						
						post4_Q1,
						paste('<b>', post4_median, '</b>', sep=''),
						post4_Q2

					),
					line = list(color = dark_grey),
					fill = list(color = c(light_grey, col_post1, col_post1, col_post1, col_post4, col_post4, col_post4)),
					align = c('left', 'center'),
					font = list(color = c(green, dark_grey,  size = 30))
				), 
				width = 0.75*as.numeric(input$dimension[1]), height = 0.75*as.numeric(input$dimension[2])
			) %>% layout(title=title_plot) %>%
				layout(plot_bgcolor='rgb(1,1,1,0)') %>% 
				layout(paper_bgcolor='rgb(1,1,1,0)') 
			
			return(table_estimations)
		}
	})
	
	## READ THE GOODNESS OF FIT TEST (GOF)
	gof_table <- reactive({
	fileName = input$results
	if (is.null(fileName)){
		return(NULL)
	}else{
		if( input$posterior_choice == 1){
			# if interested by the posterior
			x = allData()[['gof_table']]
		}else{
			if( input$posterior_choice == 2){
				# if interested by the optimized posterior
				x = allData()[['gof2_table']]
			}
		}
		return(x)
		}
	})
	
	output$gof_table = DT::renderDataTable(
		if(is.null(input$results)){
			return(NULL)
		}
		else{
			datatable(gof_table(), options = list(pageLength = 40)) %>% formatStyle('pvals_fdr_corrected', target = 'row', backgroundColor = styleInterval(cuts=c(0.01, 0.05), values=c("#fc9272", "#fee0d2", "#99d8c9")))
		}
	)
	
	# display the gof table 
	output$display_gof_table <- renderUI({
	if(is.null(input$results)){
		return(NULL)
	}else{
		DT::dataTableOutput("gof_table")
		}
	})

	
	# distribution summary statistics
	output$density_statistic <- renderPlotly({
		fileName = input$results
		if (is.null(fileName)){
			return(NULL)
		}else{
			stat_name = input$stat_name
		
			rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
			theme_set(theme_classic())
			figure = list()
			x = allData()[['distribution_PCA']]
			
			toRemove = c(1)
			for(i in 1:(ncol(x)-1)){ # -1 to keep the column prior/posterior/optimized/observed
				if(sd(x[,i]) < 0.00001){
					toRemove = c(toRemove, i)
				}
			}
			toRemove = c(toRemove, grep('max', colnames(x)))
			toRemove = c(toRemove, grep('min', colnames(x)))
			toRemove = unique(toRemove)
			x = x[, -toRemove]


			# keep simulations
			simulations = select(filter(x, origin%in%c('prior', 'posterior', 'optimized posterior')), c('origin', stat_name))
			colnames(simulations) = c('distribution', 'stat')

			# keep observation 
			observation = select(filter(x, origin%in%c('observed dataset')), c('origin', stat_name))
			colnames(observation) = c('distribution', 'stat')


			p <- ggplot(simulations, aes(stat, fill = distribution)) +
			#	geom_histogram( alpha = 0.7 ) +
				geom_density(alpha = 0.7, size = 0.25) +
				scale_fill_manual(values=c("darkgray", viridis_pal(option="D")(5)[c(1,4)])[3:1]) +
				geom_vline(xintercept=observation$stat, lwd=1.2, col=viridis_pal(option="D")(5)[5]) +
				theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), legend.text = element_text(size = 15), panel.background = element_rect(fill = '#ffffff')) +
				scale_x_continuous(name = stat_name)

			return(p)
		}
	})


	## plot the SFS : get the table
	table_sfs <- reactive({
		fileName = input$results
		if(is.null(fileName)){
			return (NULL)
		}else{
			if( input$posterior_choice == 1){
				# if interested by the posterior
				table_sfs = allData()[['gof_sfs']]
			}else{
				if( input$posterior_choice == 2){
					# if interested by the first optimized posterior
					table_sfs = allData()[['gof2_sfs']]
				}
				
			}

			if(allData()[['users_infos']][1,2]==2){
				# if nPops = 2
				## inelegant way to remove 'useless' bins...
				table_sfs[4,which(log10(table_sfs[1,])==-Inf)] = NA
				table_sfs[3,which(log10(table_sfs[1,])==-Inf)] = NA
				## end ot the inelegant block
			}
			
			return(table_sfs)
		}
	})
	
		
	## plot the SFS : display the matrix
	output$sfs_observed_2pops <- renderPlotly({
		if( is.null(input$results)){
			return(NULL)
		}else{
			f=list(
				family = "Arial",
				size = 20,
				color = "black"
			)

			f2=list(
				family = "Arial",
				size = 16,
				color = "black"
			)

			f_legend=list(
				family = "Arial",
				size = 16,
				color = "black",
				color = "#000"
			)

			xlab = list(
				title=allData()[['users_infos']][2,2],
				titlefont=f,
				tickfont=f2
			)

			ylab = list(
				title=allData()[['users_infos']][3,2],
				titlefont=f,
				tickfont=f2
			)
			
			nameA = allData()[['users_infos']][2,2]
			nameB = allData()[['users_infos']][3,2]
			noms=matrix(unlist(strsplit(names(table_sfs()), '_')), byrow=T, ncol=2)
			dat = data.frame(x=as.numeric(substr(noms[,1], 3, 10)), y=as.numeric(substr(noms[,2], 3, 10)), z=log10(as.numeric(table_sfs()[1,]))) # f(A) // f(B) // nSNPs
			dat$z[dat$z==-Inf] = NA
			dat$z[which(dat[,1]==0 & dat[,2]==1 | dat[,1]==1 & dat[,2]==0)] = NA
			
			plot_obs = plot_ly(width = 0.75*as.numeric(input$dimension[1])/2, height = 0.42*as.numeric(input$dimension[2]), colors=rev(viridis_pal(option='D')(100))) %>%
			  layout(autosize = FALSE,
				legend=list(orientation = 'h', y=1.05, font=f_legend),
				hoverlabel = list(font=list(size=20, color='#C7F464'), bordercolor='#556270', bgcolor='#556270'),
				annotations = list( text = 'observed jSFS (log10(nSNPs))', font=list(size=22), xanchor = "center", yanchor="bottom", xref="paper", yref="paper", align="center", showarrow=FALSE, x=0.5, y=1),
				autosize = T, margin = list(l=50, r=50, b=80, t=40, pad=2),
				xaxis=xlab, yaxis=ylab) %>%
			  add_trace(data = dat, x = ~x, y = ~y, z = ~z, type = "heatmap",
				    hoverinfo = 'text',
				    text = ~paste(paste(nameA, ': ', sep=''), dat$x,
						  paste('<br>', nameB, ': ', sep=''), dat$y,
						  paste('<br>number of observed SNPs: ', 10**dat$z, sep=''))
				)
			return(plot_obs)

		}
	})

	## plot the SFS : display the matrix
	output$sfs_expected_2pops <- renderPlotly({
		if( is.null(input$results)){
			return(NULL)
		}else{
			f=list(
				family = "Arial",
				size = 20,
				color = "black"
			)

			f2=list(
				family = "Arial",
				size = 16,
				color = "black"
			)

			f_legend=list(
				family = "Arial",
				size = 16,
				color = "black",
				color = "#000"
			)

			xlab = list(
				title=allData()[['users_infos']][2,2],
				titlefont=f,
				tickfont=f2
			)

			ylab = list(
				title=allData()[['users_infos']][3,2],
				titlefont=f,
				tickfont=f2
			)
			
			nameA = allData()[['users_infos']][2,2]
			nameB = allData()[['users_infos']][3,2]
			noms=matrix(unlist(strsplit(names(table_sfs()), '_')), byrow=T, ncol=2)
			dat = data.frame(x=as.numeric(substr(noms[,1], 3, 10)), y=as.numeric(substr(noms[,2], 3, 10)), z=log10(as.numeric(table_sfs()[2,]))) # f(A) // f(B) // nSNPs
			dat$z[dat$z==-Inf] = NA
			dat$z[which(dat[,1]==0 & dat[,2]==1 | dat[,1]==1 & dat[,2]==0)] = NA
			
			plot_exp = plot_ly(width = 0.75*as.numeric(input$dimension[1])/2, height = 0.42*as.numeric(input$dimension[2]), colors=rev(viridis_pal(option='D')(100))) %>%
			  layout(autosize = FALSE,
				legend=list(orientation = 'h', y=1.05, font=f_legend),
				hoverlabel = list(font=list(size=20, color='#C7F464'), bordercolor='#556270', bgcolor='#556270'),
				annotations = list( text = 'expected jSFS (log10(nSNPs))', font=list(size=22), xanchor = "center", yanchor="bottom", xref="paper", yref="paper", align="center", showarrow=FALSE, x=0.5, y=1),
				autosize = T, margin = list(l=50, r=50, b=80, t=40, pad=2),
				xaxis=xlab, yaxis=ylab) %>%
			  add_trace(data = dat, x = ~x, y = ~y, z = ~z, type = "heatmap",
				    hoverinfo = 'text',
				    text = ~paste(paste(nameA, ': ', sep=''), dat$x,
						  paste('<br>', nameB, ': ', sep=''), dat$y,
						  paste('<br>number of expected SNPs: ', 10**dat$z, sep=''))
				)
			return(plot_exp)
		}
	})
	
	## plot the SFS : display the matrix
	output$sfs_diff_2pops <- renderPlotly({
		if( is.null(input$results)){
			return(NULL)
		}else{
			f=list(
				family = "Arial",
				size = 20,
				color = "black"
			)

			f2=list(
				family = "Arial",
				size = 16,
				color = "black"
			)

			f_legend=list(
				family = "Arial",
				size = 16,
				color = "black",
				color = "#000"
			)

			xlab = list(
				title=allData()[['users_infos']][2,2],
				titlefont=f,
				tickfont=f2
			)

			ylab = list(
				title=allData()[['users_infos']][3,2],
				titlefont=f,
				tickfont=f2
			)
			
			nameA = allData()[['users_infos']][2,2]
			nameB = allData()[['users_infos']][3,2]
			noms=matrix(unlist(strsplit(names(table_sfs()), '_')), byrow=T, ncol=2)
			dat = data.frame(x=as.numeric(substr(noms[,1], 3, 10)), y=as.numeric(substr(noms[,2], 3, 10)), z=as.numeric(table_sfs()[3,])) # f(A) // f(B) // nSNPs
			dat$z[dat$z==-Inf] = NA
			dat$z[which(dat[,1]==0 & dat[,2]==1 | dat[,1]==1 & dat[,2]==0)] = NA
			
			plot_diff = plot_ly(width = 0.75*as.numeric(input$dimension[1])/2, height = 0.42*as.numeric(input$dimension[2]), colors=rev(viridis_pal(option='D')(100))) %>%
			  layout(autosize = FALSE,
				legend=list(orientation = 'h', y=1.05, font=f_legend),
				hoverlabel = list(font=list(size=20, color='#C7F464'), bordercolor='#556270', bgcolor='#556270'),
				annotations = list( text = 'expected - observed (nSNPs)', font=list(size=22), xanchor = "center", yanchor="bottom", xref="paper", yref="paper", align="center", showarrow=FALSE, x=0.5, y=1),
				autosize = T, margin = list(l=50, r=50, b=80, t=40, pad=2),
				xaxis=xlab, yaxis=ylab) %>%
			  add_trace(data = dat, x = ~x, y = ~y, z = ~z, type = "heatmap",
				    hoverinfo = 'text',
				    text = ~paste(paste(nameA, ': ', sep=''), dat$x,
						  paste('<br>', nameB, ': ', sep=''), dat$y,
						  paste('<br>exp-obs = ', dat$z, sep=''))
				)
			return(plot_diff)
		}
	})
	
	## plot the SFS : display the matrix
	output$sfs_pval_2pops <- renderPlotly({
		if( is.null(input$results)){
			return(NULL)
		}else{
			f=list(
				family = "Arial",
				size = 20,
				color = "black"
			)

			f2=list(
				family = "Arial",
				size = 16,
				color = "black"
			)

			f_legend=list(
				family = "Arial",
				size = 16,
				color = "black",
				color = "#000"
			)

			xlab = list(
				title=allData()[['users_infos']][2,2],
				titlefont=f,
				tickfont=f2
			)

			ylab = list(
				title=allData()[['users_infos']][3,2],
				titlefont=f,
				tickfont=f2
			)
			
			nameA = allData()[['users_infos']][2,2]
			nameB = allData()[['users_infos']][3,2]
			noms=matrix(unlist(strsplit(names(table_sfs()), '_')), byrow=T, ncol=2)
			dat = data.frame(x=as.numeric(substr(noms[,1], 3, 10)), y=as.numeric(substr(noms[,2], 3, 10)), z=as.numeric(table_sfs()[4,])) # f(A) // f(B) // nSNPs
			dat$z[dat$z==-Inf] = NA
			dat$z[which(dat[,1]==0 & dat[,2]==1 | dat[,1]==1 & dat[,2]==0)] = NA
			
			plot_diff = plot_ly(width = 0.75*as.numeric(input$dimension[1])/2, height = 0.42*as.numeric(input$dimension[2]), colors=rev(viridis_pal(option='D')(100)), zmin=0, zmax=0.5) %>%
			  layout(autosize = FALSE,
				legend=list(orientation = 'h', y=1.05, font=f_legend),
				hoverlabel = list(font=list(size=20, color='#C7F464'), bordercolor='#556270', bgcolor='#556270'),
				annotations = list( text = 'p-values', font=list(size=22), xanchor = "center", yanchor="bottom", xref="paper", yref="paper", align="center", showarrow=FALSE, x=0.5, y=1),
				autosize = T, margin = list(l=50, r=50, b=80, t=40, pad=2),
				xaxis=xlab, yaxis=ylab) %>%
			  add_trace(data = dat, x = ~x, y = ~y, z = ~z, type = "heatmap",
				    hoverinfo = 'text',
				    text = ~paste(paste(nameA, ': ', sep=''), dat$x,
						  paste('<br>', nameB, ': ', sep=''), dat$y,
						  paste('<br>p-value = ', dat$z, sep=''))
				)
			return(plot_diff)
		}
	})
	

	output$sfs_observed_1pop <- renderPlotly({
		if( is.null(input$results)){
			return(NULL)
		}else{
			f=list(
				family = "Arial",
				size = 26,
				color = "black"
			)

			f2=list(
				family = "Arial",
				size = 20,
				color = "black"
			)

			f_legend=list(
				family = "Arial",
				size = 20,
				color = "black",
				color = "#000"
			)

			xlab = list(
				title='Frequency',
				titlefont=f,
				tickfont=f2
			)

			ylab = list(
				title='Number of SNPs',
				titlefont=f,
				tickfont=f2
			)
			
			nClasses <- as.numeric(matrix(unlist(strsplit(colnames(table_sfs()), 'A')), ncol=2, byrow=T)[,2])
			observed <- as.numeric(table_sfs()[1,])
			expected <- as.numeric(table_sfs()[2,])
			
			data_sfs <- data.frame(nClasses, observed, expected)
			data_sfs_reshape <- data_sfs %>%
			  gather(sfs, Count, observed:expected)

			if(input$posterior_choice == 1){ # if posterior
				couleurs = viridis_pal(option = "D")(5)[c(1,5)]
			}else{ # if optimized posterior
				couleurs = viridis_pal(option = "D")(5)[c(4,5)]
			}
			data_sfs_reshape %>%
				plot_ly(type = "bar",
					x = ~nClasses,
					y = ~Count,
					color = ~sfs,
					colors = couleurs, width = (0.75*as.numeric(input$dimension[1])), height = 0.65*as.numeric(input$dimension[2])) %>%
					layout(xaxis=xlab, yaxis=ylab, legend=list(orientation = 'h', y=1.05, font=f_legend), hoverlabel = list(font=list(size=20)))
		}
	})
	
	output$parameters_estimates <- renderUI({
		if(is.null(input$results)){
			return(NULL)
		}else{
			if(allData()[['users_infos']][1,2]==2){
				# if nSpecies == 2
				tabsetPanel(
					tabPanel('Posterior', uiOutput('output_posterior_2pops')),
					tabPanel('Highest Posterior Density',
						fluidPage(
						selectInput('modeInference', label = 'prediction method', choices = list('neural network' = 'neural network', 'random forest' = 'random forest'), selected = 'neural network'),
						plotlyOutput(outputId = 'table_parameters_2pops')
						)
					),
					tabPanel("PCA", selectInput("PCA_parameters_choice", label = h4("PCA on parameters"), choices = list("Plot" = 1, "Table" = 2), selected = 1), hr(), uiOutput("display_PCA_parameter"))
				)
			}else{
				if(allData()[['users_infos']][1,2]==1){
					# if nSpecies == 2
					tabsetPanel(
						tabPanel("Posterior", uiOutput("output_posterior_1pop")),
						tabPanel("Highest Posterior Density", selectInput('modeInference', label = 'prediction method', choices = list('neural network' = 'neural network', 'random forest' = 'random forest'), selected = 'neural network'), plotlyOutput(outputId = "table_parameters_1pop")),
						tabPanel("PCA", selectInput("PCA_parameters_choice", label = h4("PCA on parameters"), choices = list("Plot" = 1, "Table" = 2), selected = 1), hr(), uiOutput("display_PCA_parameter"))
					)
				}
			}
		}
	})
	
	output$gof <- renderUI({
		if(is.null(input$results)){
			return(NULL)
		}else{
			if(allData()[['users_infos']][1,2]==2){
				# if nSpecies == 2
				tabsetPanel(
					tabPanel("Distribution of statistics", uiOutput("distribution_gof")),
					tabPanel("P-values", selectInput("posterior_choice", label = h4("Select the parameters estimate"), choices = list("Posterior" = 1, "Optimized posterior" = 2), selected = 1), hr(), uiOutput("display_gof_table")),
					tabPanel("SFS", h4(textOutput("selected_output")), uiOutput("display_sfs_table")),
					tabPanel("PCA", selectInput("PCA_gof_choice", label = h4("PCA on summary statistics"), choices = list("Plot" = 1, "Table" = 2), selected = 1), hr(), uiOutput("display_PCA_gof"))
				)
			}else{
				if(allData()[['users_infos']][1,2]==1){
					# if nSpecies == 1
					tabsetPanel(
						tabPanel("Distribution of statistics", uiOutput("distribution_gof")),
						tabPanel("P-values", selectInput("posterior_choice", label = h4("Select the parameters estimate"), choices = list("Posterior" = 1, "Optimized posterior" = 2), selected = 1), hr(), uiOutput("display_gof_table")),
						tabPanel("SFS", h4(textOutput("selected_output")), uiOutput("display_sfs_table")),
						tabPanel("PCA", selectInput("PCA_gof_choice", label = h4("PCA on summary statistics"), choices = list("Plot" = 1, "Table" = 2), selected = 1), hr(), uiOutput("display_PCA_gof"))
					)
				}

			}
		}
	})


	output$selected_output <- renderText({
		if( input$posterior_choice == 1 ){
			posterior = 'posterior'
		}else{
			if( input$posterior_choice == 2 ){
				posterior = 'optimized posterior'
			}
		}
		paste('Selected parameters estimate', posterior, sep=' : ' )
	})

	
	output$display_PCA_gof <- renderUI({
		if(is.null(input$PCA_gof_choice)){
			return(NULL)
		}else{
			if(input$PCA_gof_choice == 1){
				fluidPage(
					fluidRow( width = 12, style="margin-top:-3em",
						column(3, selectInput("axe1", label = h4("x-axis"), choices = list("PC1" = 1, "PC2" = 2, "PC3" = 3), selected = 1)),
						column(3, selectInput("axe2", label = h4("y-axis"), choices = list("PC1" = 1, "PC2" = 2, "PC3" = 3), selected = 2))
					),
					fluidRow( width = 12, plotlyOutput(outputId = "plotly_PCA_gof_2D"))
				)
			}else if(input$PCA_gof_choice == 2){
				fluidRow( width = 12, plotlyOutput(outputId = "table_PCA_gof"))
			}
		}
	})

	
	output$display_PCA_parameter <- renderUI({
		if(is.null(input$PCA_parameters_choice)){
			return(NULL)
		}else{
			if(input$PCA_parameters_choice == 1){
				fluidRow( width = 12, plotlyOutput(outputId = "plot_PCA_parameters"))
			}else if(input$PCA_parameters_choice == 2){
				fluidRow( width = 12, plotlyOutput(outputId = "table_PCA_parameters"))
			}
		}
	})

	
	output$plot_PCA_parameters <- renderPlotly({
		fileName = input$results
		
		if (is.null(fileName)){
			return(NULL)
		}else{
			green = "#C7F464"
			dark_grey = "#1e2b37"
			light_grey = "#556270"

			x = allData()[['posterior']]
			y3 = allData()[['optimized_posterior']]
		
			origin = c(rep("posterior", nrow(x)), rep("optimized posterior", nrow(y3)))

			posterior = which(origin == "posterior")
			optimized3 = which(origin == "optimized posterior")

#			data = rbind(x, y, y2, y3)
			data = rbind(x, y3)
			data = cbind(data, origin)

			res.pca <- PCA(data[, -ncol(data)], graph = FALSE, ncp=3)

			l <- list( font = list( family = "sans-serif", size = 18 ), orientation = 'v' )

			trace1 <- list(
				mode = "markers", 
				name = "posterior", 
				type = "scatter3d", 
				x = res.pca$ind$coord[,1][posterior],
				y = res.pca$ind$coord[,2][posterior],
				z = res.pca$ind$coord[,3][posterior]
			)

			trace2 <- list(
				mode = "markers", 
				name = "optimized posterior", 
				type = "scatter3d", 
				x = res.pca$ind$coord[,1][optimized3],
				y = res.pca$ind$coord[,2][optimized3],
				z = res.pca$ind$coord[,3][optimized3]
			)
			
			layout <- list(
				scene = list(
					xaxis = list(title = paste("PC1 (", round(res.pca$eig[,2][1], 2), "%)", sep=''), showline = FALSE), 
					yaxis = list(title = paste("PC2 (", round(res.pca$eig[,2][2], 2), "%)", sep=''), showline = FALSE), 
					zaxis = list(title = paste("PC3 (", round(res.pca$eig[,2][3], 2), "%)", sep=''), showline = FALSE)
				), 
				title = "PCA of parameters (3D)"
			)

			p1 <- plot_ly(type = 'scatter', mode = 'markers', width = (0.75*as.numeric(input$dimension[1])), height = 0.5*as.numeric(input$dimension[2])) %>%
				add_trace( mode=trace1$mode, name=trace1$name, type=trace1$type, x=trace1$x, y=trace1$y, z=trace1$z, marker = list(size = 10, color = viridis_pal(option='D')(5)[1]), alpha=0.8) %>%
				add_trace( mode=trace2$mode, name=trace2$name, type=trace2$type, x=trace2$x, y=trace2$y, z=trace2$z, marker = list(size = 11, color = viridis_pal(option='D')(5)[4])) %>%
				layout( scene=layout$scene, title=layout$title, legend=l, xaxis = list(showticklabels=F, zeroline=F, showline=F, showgrid=F), yaxis = list(showticklabels=F, zeroline=F, showline=F, showgrid=F), legend=list(size=12) )
			
			return( p1 )
		}
	})


	output$table_PCA_parameters <- renderPlotly({
		fileName = input$results
		
		if (is.null(fileName)){
			return(NULL)
		}else{
			green = "#C7F464"
			dark_grey = "#1e2b37"
			light_grey = "#556270"
			
			x = allData()[['posterior']]#read.table(paste(rootName, "/best_model/posterior_bestModel.txt", sep=''), h=T)
			y3 = allData()[['optimized_posterior']]#read.table(paste(rootName, "/best_model_5/posterior_bestModel.txt", sep=''), h=T)
			
			origin = c(rep("posterior", nrow(x)), rep("optimized posterior", nrow(y3)))


			data = rbind(x, y3)
			data = cbind(data, origin)

			res.pca <- PCA(data[, -ncol(data)], graph = FALSE, ncp=3)

			x = t(round(res.pca$var$contrib, 2))
			x = rbind( paste('<b>', rownames(res.pca$var$contrib), '</b>', sep=''), x)
			p1 <- plot_ly(type = 'table',
				header = list(values = c('<b>Parameters</b>', '<b>Dim. 1</b>', '<b>Dim. 2</b>', '<b>Dim. 3</b>'), line = list(color = dark_grey), fill = list(color = dark_grey), align = c('left','center'), font = list(color = green, size = 15)),
				cells = list( values=x, line = list(color = dark_grey), fill = list(color = c(light_grey, 'GhostWhite')), align = c('left','center'), font = list(color = c(green, dark_grey), size = 15)), width = (0.75*as.numeric(input$dimension[1])), height = 0.5*as.numeric(input$dimension[2])
			)
		
			return( p1 )
		}
	})

	
	output$plotly_PCA_gof_2D <- renderPlotly({
		fileName = input$results
		
		if (is.null(fileName)){
			return(NULL)
		}else{
			nspecies = allData()[['users_infos']] #read.csv(paste(rootName, "/general_infos.txt", sep=''), h=F)
			coord_PCA_SS = allData()[['coord_PCA_SS']] #read.table(paste(rootName, "/table_coord_PCA_SS.txt", sep=''), h=T, sep='\t')
			contrib_PCA_SS = allData()[['contribution_PCA']]#read.table(paste(rootName, "/table_contrib_PCA_SS.txt", sep=''), h=T, sep='\t')
			eigen = allData()[['eigen']]#read.table(paste(rootName, "/table_eigenvalues_PCA_SS.txt", sep=''), h=T, sep='\t')
				
			observed = which(coord_PCA_SS$origin == 'observed dataset')
			prior = which(coord_PCA_SS$origin == 'prior')
			posterior = which(coord_PCA_SS$origin == 'posterior')

			if(nspecies[1,2] == 2){	
				optimized_posterior = which(coord_PCA_SS$origin == 'optimized posterior')
			}else{
				optimized_posterior = which(coord_PCA_SS$origin == 'optimized posterior')
			}

			axe1 = as.numeric(input$axe1)
			axe2 = as.numeric(input$axe2)
			trace1 <- list(
				mode = "markers", 
				name = "prior", 
				type = "scatter", 
				x = coord_PCA_SS[,axe1][prior],
				y = coord_PCA_SS[,axe2][prior]
			)

			trace2 <- list(
				mode = "markers", 
				name = "posterior", 
				type = "scatter", 
				x = coord_PCA_SS[,axe1][posterior],
				y = coord_PCA_SS[,axe2][posterior]
			)

			trace3 <- list(
				mode = "markers", 
				name = "optimized posterior", 
				type = "scatter", 
				x = coord_PCA_SS[,axe1][optimized_posterior],
				y = coord_PCA_SS[,axe2][optimized_posterior]
			)

			trace6 <- list(
				mode = "markers", 
				name = "observed dataset", 
				x = coord_PCA_SS[,axe1][observed],
				y = coord_PCA_SS[,axe2][observed]
			)

			l <- list( font = list( family = "sans-serif", size = 19 ), orientation = 'v', marker = list( size = c(30,30,30,30,30,30) ))

			layout <- list(
				scene = list(
					xaxis = list(title = paste("PC",axe1, " (", round(eigen[,2][axe1], 2), "%)", sep=''), showline = T), 
					yaxis = list(title = paste("PC",axe2, " (", round(eigen[,2][axe2], 2), "%)", sep=''), showline = T)
				), 
				title = "PCA of goodness-of-fit (3D)"
			)

			xaxis = list(title = paste("PC", axe1, " (", round(eigen[,2][axe1], 2), "%)", sep=''), showline = T)
			yaxis = list(title = paste("PC", axe2, " (", round(eigen[,2][axe2], 2), "%)", sep=''), showline = T)

			#p <- plot_ly(type = 'scatter', mode = 'markers', width = (0.75*as.numeric(input$dimension[1])), height = 0.65*as.numeric(input$dimension[2])) %>%
			p <- plot_ly(type = 'scatter', mode = 'markers', width = (0.5*as.numeric(input$dimension[1])), height = 0.5*as.numeric(input$dimension[2])) %>%
				add_trace(type = 'scatter', mode=trace1$mode, name=trace1$name, type=trace1$type, x=trace1$x, y=trace1$y, marker = list(size = 8, color='darkgray')) %>%
				add_trace(type = 'scatter', mode=trace2$mode, name=trace2$name, type=trace2$type, x=trace2$x, y=trace2$y, marker = list(size = 8, color = viridis_pal(option='D')(5)[1])) %>%
				add_trace(type = 'scatter', mode=trace3$mode, name=trace3$name, type=trace3$type, x=trace3$x, y=trace3$y, marker = list(size = 8, color = viridis_pal(option='D')(5)[4])) %>%
				add_trace(type = 'scatter', mode=trace6$mode, name=trace6$name, type=trace6$type, x=trace6$x, y=trace6$y, marker = list(size = 10, color = viridis_pal(option='D')(5)[5])) %>%
				layout(xaxis = xaxis, yaxis = yaxis, font = list( family = "sans-serif", size = 19 ))

					return(p)
			}})
			
			output$table_PCA_gof <- renderPlotly({
				fileName = input$results
				
				if (is.null(fileName)){
					return(NULL)
				}else{
					green = "#C7F464"
					dark_grey = "#1e2b37"
					light_grey = "#556270"
					
					data = allData()[['contribution_PCA']]#read.table( paste(rootName, "/table_contrib_PCA_SS.txt", sep=''), h=T, sep='\t')
					
					data = data[ order(data[,1], decreasing=T), ]
						
					
					x = t(round(data, 2))
					x = rbind( paste('<b>', rownames(data), '</b>', sep=''), x)
					p1 <- plot_ly(type = 'table',
						header = list(values = c('<b>Parameters</b>', '<b>Dim. 1</b>', '<b>Dim. 2</b>', '<b>Dim. 3</b>'), line = list(color = dark_grey), fill = list(color = dark_grey), align = c('left','center'), font = list(color = green, size = 15)),
						cells = list( values=x, line = list(color = dark_grey), fill = list(color = c(light_grey, 'GhostWhite')), align = c('left','center'), font = list(color = c(green, dark_grey), size = 15)), width = (0.75*as.numeric(input$dimension[1])), height = 0.5*as.numeric(input$dimension[2])
						#cells = list( values=x, line = list(color = dark_grey), fill = list(color = c(light_grey, 'GhostWhite')), align = c('left','center'), font = list(color = c(green, dark_grey), size = 15))
					)
			return( p1 )

	}})


	
	output$display_sfs_table <- renderUI({
		if(is.null(table_sfs())){
			return(NULL)
		}else{
			if(allData()[['users_infos']][1,2]==2){
			# if nSpecies == 2
				fluidPage(
					hr(),
					
					fluidRow( style="margin-top:-4em",
						width = 12,
						column(width=6, offset = 0, style='padding:30px;',
							plotlyOutput(outputId = "sfs_observed_2pops")
						),
						column(width=6, offset = 0, style='padding:30px;',
							plotlyOutput(outputId = "sfs_expected_2pops")
						)
					),
					
					fluidRow( style="margin-top:-4em",
						width = 12,
						column(width=6, offset = 0, style='padding:30px;',
							plotlyOutput(outputId = "sfs_diff_2pops")
						),
						column(width=6, offset = 0, style='padding:30px;',
							plotlyOutput(outputId = "sfs_pval_2pops")
						)
					)
				)
			}else{
				if(allData()[['users_infos']][1,2]==1){
					# if nSpecies == 1
					fluidPage(
						fluidRow(
							width = 12,
							column(width=12, offset = 0, style='padding:30px;',
								plotlyOutput(outputId = "sfs_observed_1pop")
							)
						)
					)
				}
			}
		}
	})
	
	## READ THE TABLE WITH SUMMARY STATISTICS AND DEMOGRAPHICS INFERENCES
	locus_spe <- reactive({
		fileName = input$results
		
		if (is.null(fileName)){
			return(NULL)
		}else{
			infos_tmp = as.matrix(allData()[['users_infos']])#as.matrix(read.csv(paste(rootName, "/general_infos.txt", sep=''), h=F))
			nSpecies = infos_tmp[1,2]
			
			#read the table
			locus_spe_tmp = allData()[['locus_spe']]#read.table(locus_spe_name, h=T)
			locus_infos = allData()[['locus_infos']]#read.table(locus_infos_name, h=T)
			
			locus_spe = merge(locus_spe_tmp, locus_infos, by.x=1, by.y=1)
			
			# return the read object
			return(locus_spe)
		}
	})
	
 
	output$locus_spe_table <- DT::renderDataTable(
		if(is.null(input$results)) {return(loadingState())}
		else{
			columns = names(locus_spe())
			if (!is.null(input$selected_obs_statistics_to_display)){
				columns = input$selected_obs_statistics_to_display
			}
			return(locus_spe()[,columns,drop=FALSE])
		}
	)
	
	output$observed_columns_to_display <- renderUI({
		if(is.null(input$results)) {return()}
		else
			selectInput("selected_obs_statistics_to_display", "Select columns to display", names(locus_spe()), multiple = TRUE)
	})

	output$display_uploaded_results <- renderUI({
		if(is.null(input$results)) {return(loadingState())}
		else
			DT::dataTableOutput("locus_spe_table")
	})
	
 
	# GRAPH SITES
	output$plot_obs_stats_sites <- renderPlotly({
		# number of loci
		nLoci = nrow(locus_spe())
		
		f <- list(
		family = "Arial",
		size = 20
		)
		axis_x <- list(
		title = "",
		titlefont = f,
		tickfont = list(size = 20)
		)	
		axis_y <- list(
		title = "Proportion of sites",
		titlefont = f,
		tickfont = list(size = 20)
		)
		
		statistics_obs_sites = c(locus_spe()$sf_avg, locus_spe()$sxA_avg, locus_spe()$sxB_avg, locus_spe()$ss_avg)
		statistics_names_sites = rep(c("Sf", paste("Sx", allData()[['users_infos']][2,2], sep=' '), paste("Sx", allData()[['users_infos']][3,2], sep= ' '), "Ss"), each = nLoci)
		locus_names_sites = rep(locus_spe()$dataset, 4)
		
		data_obs_sites = data.frame(statistics_obs_sites, statistics_names_sites, locus_names_sites)
		
		if(input$show_points_stats_sites == T){
			graph_sites = plot_ly(data_obs_sites, y=~statistics_obs_sites, x=~statistics_names_sites, color=~statistics_names_sites, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_sites, "<br>", statistics_obs_sites), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}else{
			graph_sites = plot_ly(data_obs_sites, y=~statistics_obs_sites, x=~statistics_names_sites, color=~statistics_names_sites, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}
		
		return(graph_sites)
	})
	
	# GRAPH DIVERSITY
	output$plot_obs_stats_diversity <- renderPlotly({
		# number of loci
		nLoci = nrow(locus_spe())
		
		f <- list(
		family = "Arial",
		size = 20
		)
		axis_x <- list(
		title = "",
		titlefont = f,
		tickfont = list(size = 20)
		)
		axis_y <- list(
		title = "Index of diversity per site",
		titlefont = f,
		tickfont = list(size = 20)
		)
		statistics_obs_diversity = c(locus_spe()$piA_avg, locus_spe()$piB_avg, locus_spe()$thetaA_avg, locus_spe()$thetaB_avg)
		statistics_names_diversity = rep(c(paste("pi", allData()[['users_infos']][2,2], sep=' '), paste("pi", allData()[['users_infos']][3,2], sep=' '), paste("Watterson's theta", allData()[['users_infos']][2,2], sep= ' '), paste("Watterson's theta", allData()[['users_infos']][3,2], sep= ' ')), each = nLoci)
		locus_names_diversity = rep(locus_spe()$dataset, 4)
		
		data_obs_diversity = data.frame(statistics_obs_diversity, statistics_names_diversity, locus_names_diversity)
		
		if(input$show_points_diversity == T){
			graph_diversity = plot_ly(data_obs_diversity, y=~statistics_obs_diversity, x=~statistics_names_diversity, color=~statistics_names_diversity, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_diversity, "<br>", statistics_obs_diversity), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}else{
			graph_diversity = plot_ly(data_obs_diversity, y=~statistics_obs_diversity, x=~statistics_names_diversity, color=~statistics_names_diversity, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(4)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}
		
		return(graph_diversity)
	})
	
	
	# GRAPH TAJIMA
	output$plot_obs_stats_tajima <- renderPlotly({
		# number of loci
		nLoci = nrow(locus_spe())
		f <- list(
			family = "Arial",
			size = 20
		)
		axis_x <- list(
			title = "",
			titlefont = f,
			tickfont = list(size = 20)
		)
		axis_y <- list(
			title = "Tajima's D",
			titlefont = f,
			tickfont = list(size = 20)
		)
		statistics_obs_tajima = c(locus_spe()$DtajA_avg, locus_spe()$DtajB_avg)
		statistics_names_tajima = rep(c(paste("Tajima's D", allData()[['users_infos']][2,2], sep=' '), paste("Tajima's D", allData()[['users_infos']][3,2], sep=' ')), each = nLoci)
		locus_names_tajima = rep(locus_spe()$dataset, 2)
		data_obs_tajima = data.frame(statistics_obs_tajima, statistics_names_tajima, locus_names_tajima)
		
		if(input$show_points_Tajima==T){
			graph_tajima = plot_ly(data_obs_tajima, y=~statistics_obs_tajima, x=~statistics_names_tajima, color=~statistics_names_tajima, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_tajima, "<br>", statistics_obs_tajima), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(2)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}else{
			graph_tajima = plot_ly(data_obs_tajima, y=~statistics_obs_tajima, x=~statistics_names_tajima, color=~statistics_names_tajima, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(2)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}
		
		return(graph_tajima)
	})
	
	
	# GRAPH DIVERGENCE
	output$plot_obs_stats_divergence <- renderPlotly({
		# number of loci
		nLoci = nrow(locus_spe())
		f <- list(
			family = "Arial",
			size = 30
		)
		
		axis_x <- list(
			title = "",
			tickfont = list(size = 20)
		)
		
		axis_y <- list(
			title = "Measure of divergence/differentiation",
			titlefont = f,
			tickfont = list(size = 20)
		)
		
		statistics_obs_divergence = c(locus_spe()$divAB_avg, locus_spe()$netdivAB_avg, locus_spe()$FST_avg)
		statistics_names_divergence = rep(c("raw divergence", "net divergence", "Fst"), each = nLoci)
		locus_names_divergence = rep(locus_spe()$dataset, 3)
		data_obs_divergence = data.frame(statistics_obs_divergence, statistics_names_divergence, locus_names_divergence)
		
		if(input$show_points_divergence==T){
			graph_divergence = plot_ly(data_obs_divergence, y=~statistics_obs_divergence, x=~statistics_names_divergence, color=~statistics_names_divergence, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_divergence, "<br>", statistics_obs_divergence), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(3)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}else{
			graph_divergence = plot_ly(data_obs_divergence, y=~statistics_obs_divergence, x=~statistics_names_divergence, color=~statistics_names_divergence, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(3)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}
		
		return(graph_divergence)
	})

	
	# GRAPHE ONE POP
	output$plot_obs_stats_1pop <- renderPlotly({
		# number of loci
		nLoci = nrow(locus_spe())
		f <- list(
			family = "Arial",
			size = 30
		)
		
		axis_x <- list(
			title = "",
			tickfont = list(size = 20)
		)
		
		axis_y <- list(
			title = "Pattern of polymorphism",
			titlefont = f,
			tickfont = list(size = 20)
		)
	
		if(input$show_points_stats_sites_1pop == TRUE){	
			statistics_obs_polyM = c(locus_spe()$piA_avg, locus_spe()$thetaA_avg) # piA_avg thetaA_avg DtajA_avg
			statistics_names_polyM = rep(c(paste('pi ', as.character(allData()[['users_infos']][2,2]), sep=''), paste('theta ', as.character(allData()[['users_infos']][2,2]), sep='')), each = nLoci)
			locus_names_polyM = rep(locus_spe()$dataset, 2)
			data_obs_polyM = data.frame(statistics_obs_polyM, statistics_names_polyM, locus_names_polyM)
			graph_polyM = plot_ly(data_obs_polyM, y=~statistics_obs_polyM, x=~statistics_names_polyM, color=~statistics_names_polyM, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_polyM, "<br>", statistics_obs_polyM), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(2)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
			
			statistics_obs_TajD = c(locus_spe()$DtajA_avg) # piA_avg thetaA_avg DtajA_avg
			statistics_names_TajD = rep(paste("Tajima's D ", as.character(allData()[['users_infos']][2,2]), sep=''), each = nLoci)
			locus_names_TajD = locus_spe()$dataset
			data_obs_TajD = data.frame(statistics_obs_TajD, statistics_names_TajD, locus_names_TajD)
			graph_TajD = plot_ly(data_obs_TajD, y=~statistics_obs_TajD, x=~statistics_names_TajD, color=~statistics_names_TajD, type="box", boxpoints="all", boxmean=T, text = ~paste0("locus: ", locus_names_TajD, "<br>", statistics_obs_TajD), hoverinfo="text", width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "C")(1)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}else{
			statistics_obs_polyM = c(locus_spe()$piA_avg, locus_spe()$thetaA_avg) # piA_avg thetaA_avg DtajA_avg
			statistics_names_polyM = rep(c(paste('pi ', as.character(allData()[['users_infos']][2,2]), sep=''), paste('theta ', as.character(allData()[['users_infos']][2,2]), sep='')), each = nLoci)
			locus_names_polyM = rep(locus_spe()$dataset, 2)
			data_obs_polyM = data.frame(statistics_obs_polyM, statistics_names_polyM, locus_names_polyM)
			graph_polyM = plot_ly(data_obs_polyM, y=~statistics_obs_polyM, x=~statistics_names_polyM, color=~statistics_names_polyM, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "D")(2)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
			
			statistics_obs_TajD = c(locus_spe()$DtajA_avg) # piA_avg thetaA_avg DtajA_avg
			statistics_names_TajD = rep(paste("Tajima's D ", as.character(allData()[['users_infos']][2,2]), sep=''), each = nLoci)
			locus_names_TajD = locus_spe()$dataset
			data_obs_TajD = data.frame(statistics_obs_TajD, statistics_names_TajD, locus_names_TajD)
			graph_TajD = plot_ly(data_obs_TajD, y=~statistics_obs_TajD, x=~statistics_names_TajD, color=~statistics_names_TajD, type="box", boxpoints="none", boxmean=T, width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2]), colors = viridis_pal(option = "C")(1)) %>% layout(xaxis = axis_x, yaxis = axis_y, legend=list(orientation = 'h', y=1.05, font = list(size = 15)), hoverlabel = list(font=list(size=20)) )
		}
		p <- subplot(graph_polyM, graph_TajD)
		return(p)
	})
	
	
	# LOCUS SPECIFIC MODEL COMPARISON
	output$locus_specific_model_comparison <- renderPlotly({
		f=list(
			family = "Arial",
			size = 26,
			color = "black"
		)
		
		f2=list(
			family = "Arial",
			size = 20,
			color = "black"
		)
		
		f3=list(
			family = "Arial",
			size = 14,
			color = "grey17"
		)
		
		f_legend=list(
			family = "Arial",
			size = 20,
			color = "black",
			color = "#000"
		)
		
		xlab = list(title='FST',
			titlefont=f,
			tickfont=f2
		)
		
		ylab_divergence = list(title='net divergence',
			titlefont=f,
			tickfont=f2
		)

		ylab_diversity = list(title='average pi',
			titlefont=f,
			tickfont=f2
		)
		
		lab_piA = list(title=paste('pi', allData()[['users_infos']][2,2], sep=' '),
			titlefont=f,
			tickfont=f2
		)
		
		lab_piB = list(title=paste('pi', allData()[['users_infos']][3,2], sep=' '),
			titlefont=f,
			tickfont=f2
		)
		
		threshold = input$threshold_locus_specific_model_comp
		
		allocation = as.vector(locus_spe()$allocation)
		allocation[which(locus_spe()$post_proba<threshold)] = 'ambiguous'
		y = data.frame(netdivAB=locus_spe()$netdivAB_avg, pi=(locus_spe()$piA_avg+locus_spe()$piB_avg)/2, FST=locus_spe()$FST_avg, allocation=allocation, post_prob=locus_spe()$post_prob, dataset=locus_spe()$dataset, piA=locus_spe()$piA_avg, piB=locus_spe()$piB_avg)
		
		plot_locus_modComp_2species_divergence <- y %>% plot_ly(x =~FST, y =~netdivAB, mode = 'markers', type = 'scatter', color =~allocation, legendgroup = ~allocation, colors = viridis_pal(option = "D")(3), marker= list(size=18, opacity=0.75), text = ~paste("Locus: ", dataset, '<br>Status: ', allocation, '<br>Posterior probability: ', round(post_prob, 5), '<br>Fst: ', round(FST, 5), '<br>net divergence: ', round(netdivAB, 5), paste('<br>pi', as.character(allData()[['users_infos']][2,2]), ':', sep=' '), round(piA, 5), paste('<br>pi', as.character(allData()[['users_infos']][3,2]), ':', sep=' '), round(piB, 5)),
		hoverinfo='text', width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2])) %>% layout(xaxis = xlab, yaxis = ylab_divergence, legend=list(orientation = 'h', y=1.05, font = list(size = 25), hoverlabel = list(font=list(size=20))))

		plot_locus_modComp_2species_pi_AB <- y %>% plot_ly(x =~FST, y =~pi, mode = 'markers', type = 'scatter', color =~allocation, legendgroup = ~allocation, colors = viridis_pal(option = "D")(3), showlegend=FALSE, marker= list(size=18, opacity=0.75), text = ~paste("Locus: ", dataset, '<br>Status: ', allocation, '<br>Posterior probability: ', round(post_prob, 5), '<br>Fst: ', round(FST, 5), '<br>net divergence: ', round(netdivAB, 5), paste('<br>pi', as.character(allData()[['users_infos']][2,2]), ':', sep=' '), round(piA, 5), paste('<br>pi', as.character(allData()[['users_infos']][3,2]), ':', sep=' '), round(piB, 5)),
		hoverinfo='text') %>% layout(xaxis = xlab, yaxis = ylab_diversity)

		plot_locus_modComp_2species_piA_piB <- y %>% plot_ly(x =~piA, y =~piB, mode = 'markers', type = 'scatter', color =~allocation, legendgroup = ~allocation, colors = viridis_pal(option = "D")(3), showlegend=FALSE, marker= list(size=12, opacity=0.65), text = ~paste("Locus: ", dataset, '<br>Status: ', allocation, '<br>Posterior probability: ', round(post_prob, 5), '<br>Fst: ', round(FST, 5), '<br>net divergence: ', round(netdivAB, 5), paste('<br>pi', as.character(allData()[['users_infos']][2,2]), ':', sep=' '), round(piA, 5), paste('<br>pi', as.character(allData()[['users_infos']][3,2]), ':', sep=' '), round(piB, 5)),
		hoverinfo='text') %>% layout(xaxis = lab_piA, yaxis = lab_piB)
		
		barplot_locus_modComp_2species <- y %>% plot_ly(x = ~allocation, color = ~allocation, mode = 'markers', type = 'scatter', legendgroup = ~allocation, colors = viridis_pal(option = "D")(3), showlegend=FALSE) %>% layout(hoverlabel = list(font=list(size=20)), yaxis= list(titlefont=f, tickfont=f2), xaxis = list(titlefont=f, tickfont=f2))
		# width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2])

		figure = subplot( plot_locus_modComp_2species_divergence, barplot_locus_modComp_2species, plot_locus_modComp_2species_pi_AB, plot_locus_modComp_2species_piA_piB, nrows = 2, widths=c(3/4, 1/4), shareX = FALSE, titleX = TRUE, titleY = TRUE)
		
		return(figure)
	})


	output$page_greyzone <- renderUI({
		fileName = input$results
			if (is.null(fileName)){
				htmltools::div(style = "display:inline-block", plotlyOutput("plot_greyzone", width = "auto"))
			}else{
				rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
				fluidPage(style="margin-top:-3em",
					if( rootName%in%allData()[['meta']][,1]==FALSE ){
						fluidRow(
							HTML('<H4>Clicking on this button <b>will save</b>:<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>1)</b> the user&#39;s email address to contact him/her for future collaborative meta-analysis<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>2)</b> the names of the organisms and the position of the point on the graph<br><b>Everything else uploaded by the user <u>will be deleted</u> from the server</b>.<br><br>An unfortunate click <b>can be cancelled</b> by uploading the same archive a second time and then clicking on <b>REMOVE THE POINT</b> button</H4>'),
							actionButton("update_greyzone", "UPDATE THE FIGURE WITH YOUR RESULTS")
						)
					}else{
						fluidRow(
							HTML('<H4><b>This analysis is already part of the figure</b>.<br>You can remove it by clicking on the <b>REMOVE THE POINT</b> button.</H4>'),
							actionButton("downgrade_greyzone", "REMOVE THE POINT")
						)

					},
					
					fluidRow(
						htmltools::div(style = "display:inline-block", plotlyOutput("plot_greyzone", width = "auto"))
					)
				)
			}
	})
	
	observeEvent(input$update_greyzone, {if (input$update_greyzone == 1)	removeUI(selector='#update_greyzone', immediate=TRUE)}, autoDestroy=TRUE)
	
	observeEvent(input$downgrade_greyzone, {if (input$downgrade_greyzone == 1)	removeUI(selector='#downgrade_greyzone', immediate=TRUE)}, autoDestroy=TRUE)
	
	observeEvent(input$update_greyzone,
		{if (input$update_greyzone == 1)
			fileName = input$results
			if (is.null(fileName)){
			}else{	
				rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
				species_A_user = as.character(allData()[['users_infos']][2,2])
				species_B_user = as.character(allData()[['users_infos']][3,2])
				author_user = as.character(allData()[['users_infos']][5,2])
				modelComp = allData()[['hierarchical']]
				ABCstat = allData()[['ABCstatGlobal']] 
			
				best_model = as.character(modelComp[2,1])
				probaMigration = as.numeric(as.character(modelComp[3,1]))
				if(best_model == 'isolation'){
					probaMigration = 1 - probaMigration
					pMigHetero = 0
					status = 'species'
				}else{
					pMigHetero = as.numeric(as.character(modelComp[3,3]))
					seuil1 = 0.6419199
					if(pMigHetero >= seuil1){
						status = 'semi-isolated species'
					}else{
						status = 'populations'
					}
				}
				
				piA_user = ABCstat$piA_avg
				piB_user = ABCstat$piB_avg
				divergence_user = ABCstat$netdivAB_avg

				res = read.table('metaanalysis.txt', sep='\t', h=T)
				
				if(sum(rootName%in%res$yaml) == 0){ # if the rootName is absent from the metaanalysis
					
					colnames = c('yaml', 'mail_address', 'speciesA', 'speciesB', 'bestModel', 'probaMigration', 'probaMigHetero', 'status', 'piA', 'piB', 'netDivergence')
					a = c(rootName, author_user, species_A_user, species_B_user, best_model, probaMigration, pMigHetero, status, piA_user, piB_user, divergence_user)
			
					obs = data.frame(yaml=a[1], mail_address=a[2], speciesA=a[3], speciesB=a[4], bestModel=a[5], probaMigration=a[6], probaMigHetero=a[7], status=a[8], piA=a[9], piB=a[10], netDivergence=a[11])
					res = rbind(res, obs)
					write.table(res, 'metaanalysis.txt', col.names=colnames, row.names=F, sep='\t', quote=F)
				}
			}
		}
	)
	
	observeEvent(input$downgrade_greyzone,
		{if (input$downgrade_greyzone == 1)
			fileName = input$results
			if (is.null(fileName)){
			}else{	
				rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
				res = read.table('metaanalysis.txt', sep='\t', h=T)
				toRemove = which(res[,1] == rootName)
				
				colnames = c('yaml', 'mail_address', 'speciesA', 'speciesB', 'bestModel', 'probaMigration', 'probaMigHetero', 'status', 'piA', 'piB', 'netDivergence')
				write.table(res[-toRemove,], 'metaanalysis.txt', col.names=colnames, row.names=F, sep='\t', quote=F)
			}
		}
	)
	
	output$plot_greyzone <- renderPlotly({
		# GREYZONE PART
		# metaanalayse
		meta = read.table('metaanalysis.txt', sep='\t', h=T)
		
		# popphyl
		x = read.table("popPhyl.txt", h=T)
		pmig_HH = x$Pongoing_migration_Mhetero_Nhetero 
		proba_migration = pmig_HH
		seuil1 = 0.6419199
		seuil2 = 0.1304469
		
		model = rep('ambiguous', length(proba_migration))
		model[which(proba_migration>=seuil1)] = "migration"
		model[which(proba_migration<seuil2)] = "isolation"
	
		divergence = x$netdivAB_avg
		divergence = log10(divergence)
	
		piA = round(x$piA_avg, 5)
		
		piB = round(x$piB_avg, 5)
		
		pattern=c("Mhetero_Nhetero", "Hetero")
		selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(x))))
		status = rep('ambiguous', nrow(x))
		heteroM = apply(x[, selectedCol], FUN="sum", MARGIN=1)
		status[which(pmig_HH>= seuil1 & heteroM >= seuil1)] = "semi-isolated species"
		status[which(pmig_HH>= seuil1 & heteroM < seuil1)] = "populations"
		status[which(pmig_HH<= seuil2)] = "species"

			
		species_A = as.character(x$spA)
		species_B = as.character(x$spB)
		
		author = rep('camille.roux@univ-lille.fr', length(species_A))
		
		if(nrow(meta) != 0){
			model = c(model, as.character(meta$bestModel))
			divergence = c(divergence, log10(as.numeric(as.character(meta$netDivergence))))
			proba_migration = c(proba_migration, as.numeric(as.character(meta$probaMigration)))
			piA = c(piA, as.numeric(as.character(meta$piA)))
			piB = c(piB, as.numeric(as.character(meta$piB)))
			status_tmp = rep('ambiguous', length(meta$status))
			status_tmp[which(meta$probaMigration>= seuil1 & meta$probaMigHetero >= seuil1)] = "semi-isolated species"
			status_tmp[which(meta$probaMigration>= seuil1 & meta$probaMigHetero < seuil1)] = "populations"
			status_tmp[which(meta$probaMigration<= seuil2)] = "species"
#			status_tmp[which(meta$probaMigration<seuil1 & meta$probaMigration>seuil2)] = as.factor('ambiguous')
			status = c(status, status_tmp)
			species_A = c(species_A, as.character(meta$speciesA))
			species_B = c(species_B, as.character(meta$speciesB))
			author = c(author, as.character(meta$mail_address))
		}
		
		# USER'S PART
		fileName = input$results
		if (is.null(fileName)){
			col = c(grey(0.25), 'turquoise', 'purple', 'red')
		}else{	
			if(allData()[['users_infos']][1,2]==2){
				col = c(grey(0.25), 'turquoise', 'purple', 'red', '#fdae61')
				
				rootName = strsplit(fileName$name, '.', fixed=T)[[1]][1]
				ABCstat = allData()[['ABCstatGlobal']]
				modelComp = allData()[['hierarchical']][-1,]
				
				divergence_user = log10(ABCstat$netdivAB_avg)
				
				model_user = as.character(modelComp[1,1])
				P = as.numeric(as.character(modelComp[2,1]))
				
				status_user = "user's point"
				
				if(model_user == 'migration'){
					proba_migration_user = P
				}else{
					proba_migration_user = 1-P
				}
				species_A_user = as.character(allData()[['users_infos']][2,2])
				species_B_user = as.character(allData()[['users_infos']][3,2])
				piA_user = ABCstat$piA_avg
				piB_user = ABCstat$piB_avg
				author_user = as.character(allData()[['users_infos']][5,2])
			
				divergence = c(divergence, divergence_user)
				model = c(model, model_user)
				status = c(status, "user's point")
				proba_migration = c(proba_migration, proba_migration_user)
				species_A = c(as.character(species_A), species_A_user)
				species_B = c(as.character(species_B), species_B_user)
				piA = c(piA, piA_user)
				piB = c(piB, piB_user)
				author = c(author, author_user)
				
			}else{
				col = c(grey(0.25), 'turquoise', 'purple', 'red')
			}
		}
		res = data.frame(divergence, model, status, proba_migration, species_A, species_B, piA, piB, author)
		
		f=list(
			family = "Arial",
			size = 26,
			color = "black"
		)
		
		f2=list(
			family = "Arial",
			size = 20,
			color = "black"
		)
		
		f_legend=list(
			family = "Arial",
			size = 20,
			color = "black",
			color = "#000"
		)
		
		xlab = list(
			title='divergence (log10)',
			titlefont=f,
			tickfont=f2,
			tickvals=c(0, -1, -2, -3, -4, -5),
			ticktext=c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
		)
		
		ylab = list(
			title='probability of ongoing migration',
			titlefont=f,
			tickfont=f2
		)
		
		p=plot_ly(data=res, x=~divergence, y=~proba_migration, type='scatter', mode = 'markers', color=~status, colors=col, marker=list(size=20), text = ~paste("species A: ", species_A, '<br>species B: ', species_B, "<br>net neutral divergence: ", round(10**divergence, 5), "<br>Probability of migration: ", round(proba_migration, 3), '<br>piA: ', piA, '<br>piB: ', piB, '<br><br>author: ', author), hoverinfo='text', width = (0.75*as.numeric(input$dimension[1])), height = 0.75*as.numeric(input$dimension[2])) %>% layout(xaxis=xlab, yaxis=ylab, legend=list(orientation = 'h', y=1.05, font=f_legend), hoverlabel = list(font=list(size=20)))
		#htmlwidgets::saveWidget(p, "figure_greyzone.html") # HTML
		#webshot::webshot("figure_greyzone.html", "figure_greyzone.pdf") # PDF -> Margin problem (cut!)
	return(p)
	})
	
	## INFORMATIONS
	### LOGOS
	output$logos <-
	renderText({
		c('<img src=https://raw.githubusercontent.com/popgenomics/ABConline/master/webinterface/pictures_folder/logos.png align="middle" height="auto" width="150%" margin="0 auto">')
	})
}

# Run the application 
options('shiny.port'=port, shiny.host=host)
shinyApp(ui = ui, server = server)

