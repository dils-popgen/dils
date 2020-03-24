# plot.ly
## requires: sudo apt-get install pandoc
library(plotly)
library(viridis)
x = read.table("popPhyl.txt", h=T)
col = c(grey(0.25), 'turquoise', 'purple', 'red')
#col = c(grey(0.75), rev(viridis(5))[3:5])
pmig_HH = x$Pongoing_migration_Mhetero_Nhetero 
proba_migration = pmig_HH
seuil1 = 0.6419199
seuil2 = 0.1304469

model = rep('ambiguous', nrow(x))
model[which(x$Pongoing_migration_Mhetero_Nhetero>=seuil1)] = "migration"
model[which(x$Pongoing_migration_Mhetero_Nhetero<seuil2)] = "isolation"

divergence = log10(x$netdivAB_avg)

piA = round(x$piA_avg, 5)
piB = round(x$piB_avg, 5)

pattern=c("Mhetero_Nhetero", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(x))))
status = rep('ambiguous', nrow(x))
heteroM = apply(x[, selectedCol], FUN="sum", MARGIN=1)
status[which(pmig_HH>= seuil1 & heteroM >= seuil1)] = "semi-isolated species"
status[which(pmig_HH>= seuil1 & heteroM < seuil1)] = "populations"
status[which(pmig_HH<= seuil2)] = "species"

species_A = x$spA
species_B = x$spB

author = rep('camille.roux@univ-lille.fr', length(species_A))
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

xlab = list(title='divergence (log10)',
	titlefont=f,
	tickfont=f2,
	tickvals=c(0, -1, -2, -3, -4, -5),
	#ticktext=c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5")))
	ticktext=c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
)

ylab = list(title='probability of ongoing migration',
	titlefont=f,
	tickfont=f2
)

p=plot_ly(data=res, x=~divergence, y=~proba_migration, color=~status, colors=col, marker=list(size=20), text = ~paste("species A: ", species_A, '<br>species B: ', species_B, "<br>net neutral divergence: ", round(10**divergence, 5), "<br>Probability of migration: ", round(proba_migration, 3), '<br>piA: ', piA, '<br>piB: ', piB, '<br><br>author: ', author), hoverinfo='text') %>% layout(xaxis=xlab, yaxis=ylab, legend=list(orientation = 'h', y=1.05, font=f_legend))
htmlwidgets::saveWidget(p, "figure_greyzone.html") # HTML
webshot::webshot("figure_greyzone.html", "figure_greyzone.pdf") # PDF -> Margin problem (cut!)

