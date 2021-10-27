
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(pscl)
library(grid)
library(cowplot)

data <- read.table("/home/fantin/z.final/9.1.DegenerationTempo/dS_FOP.txt", header = F)
colnames(data) <- c("gene", "dS", "genomicCompartment", "species", "optimal", "non_optimal")
data$FOP <- data$optimal/data$non_optimal

####For a1 mating-type haploid genomes - plot
sub_data1 <- data[c(grep("1$", data$species), grep("D$", data$species)),]

#Table with mean dS
new_data1_dS <- sub_data1 %>% 
	group_by(genomicCompartment, species) %>% 
	summarise(pip = mean(dS) ) %>% 
	spread(species, pip)

new_data1_dS <- melt(as.data.frame(new_data1_dS))
colnames(new_data1_dS) <- c("genomicCompartment", "species", "mean_dS")


#Plot smooth following logit scale
tableFOPdegeneration <- merge(sub_data1, new_data1_dS, by = c("genomicCompartment", "species"))
trialFOPdegeneration <- cbind(tableFOPdegeneration$optimal, tableFOPdegeneration$non_optimal)

modelFOPdegeneration <- glm(trialFOPdegeneration ~ mean_dS + species, 
	family = binomial,
	data = tableFOPdegeneration)


tableFOPdegeneration %>% 
	mutate(pred_logit = predict(modelFOPdegeneration)) -> tableFOPdegeneration

plotReglogDegeneration_logitScale <- tableFOPdegeneration %>% 
	ggplot() +
	aes(x = mean_dS, y = pred_logit, color = species) +
	geom_point() +
	geom_line()

pdf("plotReglogDegeneration_logitScale.pdf");
	plotReglogDegeneration_logitScale
dev.off()


#Plot smooth following natural scale
tableFOPdegeneration <- merge(sub_data1, new_data1_dS, by = c("genomicCompartment", "species"))
trialFOPdegeneration <- cbind(tableFOPdegeneration$optimal, tableFOPdegeneration$non_optimal)

modelFOPdegeneration <- glm(trialFOPdegeneration ~ mean_dS + species, 
	family = binomial,
	data = tableFOPdegeneration)

tableFOPdegeneration %>% 
	mutate(pred_prob = predict(modelFOPdegeneration, 
		type = "response")) -> tableFOPdegeneration

#plotReglogDegeneration_naturalScale <- tableFOPdegeneration %>% 
#	ggplot() +
#	geom_smooth(method = "glm",
#		aes(x = mean_dS, y = pred_prob, colour = species),
#		method.args = list(family = "binomial"),
#		se = F)


#pdf("plotReglogDegeneration_naturalScale.pdf");
#	plotReglogDegeneration_naturalScale
#dev.off()



#Table with mean FOP
new_data1_FOP <- sub_data1 %>% 
	group_by(genomicCompartment, species) %>% 
	summarise(pip = mean(FOP) ) %>% 
	spread(species, pip)

new_data1_FOP <- melt(as.data.frame(new_data1_FOP))
colnames(new_data1_FOP) <- c("genomicCompartment", "species", "mean_FOP")

#Merging...
toPLOT1 <- merge(new_data1_dS, new_data1_FOP, 
	by = c("genomicCompartment", "species") )


#Plotting...
toPLOT1$species <- revalue(toPLOT1$species,
	c("MinSalD" = "A",
	"MscKan1" = "B",
	"MviSic1" = "C",
	"MviSta1" = "D",
	"MlaSvu1" = "E",
	"MsaSof1" = "F",
	"MviSpa1" = "G",
	"MsaSac1" = "H",
	"MviSco1" = "I",
	"MviLyf1" = "J",
	"MviSin1" = "K",
	"MsdSdi1" = "L",
	"MldSil1" = "M"
	)
)

colors_species <- c(
	"A" = "lightpink",
	"B" = "lightcyan3",
	"C" = "palevioletred",
	"D"	= "chocolate",
	"E" = "lightsalmon2",
	"F" = "indianred2",
	"G" = "lightgoldenrod",
	"H" = "darkolivegreen3",
	"I" = "royalblue4",
	"J" = "gold1",
	"K" = "seagreen",
	"L" = "darkorchid2",
	"M" = "lightseagreen"
)

fill_genomicCompartment <- c(
	"AUT" = "grey80",
	"PAR" = "grey60",
	"PR_RR" = "grey50",
	"HD_RR" = "grey40",
	"purpleToCEN" = "mediumorchid3",
	"blueToCEN" = "blue3",
	"purple" = "mediumorchid1",
	"blue" = "blue1",
	"NRR" = "black",
	"orange" = "orange",
	"red" = "red",
	"pink" = "pink",
	"lightblue" = "lightblue",
	"white" = "white"
)

font_genomicCompartment <- c(
        "AUT" = "black",
        "PAR" = "black",
        "PR_RR" = "black",
        "HD_RR" = "black",
        "purpleToCEN" = "black",
        "blueToCEN" = "black",
        "purple" = "black",
        "blue" = "black",
        "NRR" = "white",
        "orange" = "black",
        "red" = "black",
        "pink" = "black",
        "lightblue" = "black",
        "white" = "black"
)


colors_all <- c(
	"A" = "lightpink",
	"B" = "lightcyan3",
	"C" = "palevioletred",
	"D" = "chocolate",
	"E" = "lightsalmon2",
	"F" = "indianred2",
	"G" = "lightgoldenrod",
	"H" = "darkolivegreen3",
	"I" = "royalblue4",
	"J" = "gold1",
	"K" = "seagreen",
	"L" = "darkorchid2",
	"M" = "lightseagreen",
	"AUT" = "black",
	"PAR" = "black",
	"NRR" = "white",
	"blueGenes" = "black",
	"purpleGenes" = "black",
	"orangeGenes" = "black",
	"lightblueStratum" = "black",
	"pinkStratum" = "black",
	"whiteStratum" = "black",
	"redGenes" = "black"
)

##########################
##########################
##########################

toPLOT1$genomicCompartment <- ordered(toPLOT1$genomicCompartment,
	levels = c("AUT", "PR_RR", "HD_RR", "PAR", "NRR", "purple", "purpleToCEN", "blue", "blueToCEN", "orange", "white", "lightblue", "pink", "red")
	)

toPLOT1$species <- ordered(toPLOT1$species, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"))

final_plot <- ggplot() +
        geom_label( aes( x = toPLOT1$mean_dS, y = toPLOT1$mean_FOP,
                label = toPLOT1$species, fill = toPLOT1$genomicCompartment, colour = toPLOT1$genomicCompartment, size = factor(toPLOT1$species) ) ) +

	scale_fill_manual(values = fill_genomicCompartment, name = "Genomic compartment\n",
	labels = c("Autosomes", "Recombining region on PR mating-type chromosome", "Recombining region on HD mating-type chromosome", "PAR",
		"Black stratum", "Purple stratum", "Non-recombining region linking the PR to the centromere", "Blue stratum",
		"Non-recombining region linking the HD to the centromere", "Orange stratum", "White stratum", "Lightblue stratum", "Pink stratum",
		"Red stratum") ) +
	scale_colour_manual(values = font_genomicCompartment) +

	scale_size_manual(values = c(3,3,3,3,3,3,3,3,3,3,3,3,3), name = "Microbotryum species\n",
	labels = c("A = M. intermedium",
        "B = M. scabiosae",
        "C = M. v. caroliniana",
        "D = M. v. tatarinowii",
        "E = M. lagerheimii",
        "F = M. saponariae",
        "G = M. v. paradoxa",
        "H = M. silenes-acaulis",
        "I = M. v. melanantha",
        "J = M. coronariae",
        "K = M. violaceum (s.s.)",
        "L = M. silenes-dioicae",
        "M = M. lychnidis-dioicae"
        )
		) +

	guides(colour = F) +
	guides(fill=guide_legend(nrow=5, title.position="top", title.hjust = 0.5, override.aes = aes(label = "")),
		size = guide_legend(nrow=3, title.position="top", title.hjust = 0.5), override.aes = aes(label = "")) +
	theme_bw() +
	xlab("Mean synonymous divergence values") +
	ylab("Mean frequency of optimal codons") +
	theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
	legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"), 
	legend.title=element_text(size = 8, face = "bold") )


legend_one <- plot_grid(get_legend(final_plot))
final_plot <- final_plot + theme(legend.position='none')

svg(paste0("/home/fantin/z.final/z.4MSfigures/Figure.4.svg"))

grid.newpage()
# Créer la mise en page : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(3, 1)))
# Une fonction pour definir une region dans la mise en page
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
# Arranger les graphiques
print(final_plot, vp = define_region(1:2, 1))
print(legend_one, vp=define_region(3, 1))


dev.off()






############################


#p1_withSE <- ggplot() +
#	geom_label(aes(x = toPLOT1$mean_dS, y = toPLOT1$mean_FOP, 
#		label = toPLOT1$species, fill = toPLOT1$genomicCompartment), size = 3) +
#	geom_smooth(method = "glm",
#		aes(x = tableFOPdegeneration$mean_dS, y = tableFOPdegeneration$pred_prob,
#		 colour = tableFOPdegeneration$species),
#		method.args = list(family = "binomial"),
#		se = T) +
#	scale_fill_manual(values = fill_genomicCompartment)
#
#p1_withoutSE <- ggplot() +
#	geom_label(aes(x = toPLOT1$mean_dS, y = toPLOT1$mean_FOP, 
#		label = toPLOT1$species, fill = toPLOT1$genomicCompartment), size = 3) +
#	geom_smooth(method = "glm",
#		aes(x = tableFOPdegeneration$mean_dS, y = tableFOPdegeneration$pred_prob,
#		 colour = tableFOPdegeneration$species),
#		method.args = list(family = "binomial"),
#		se = F) +
#	scale_fill_manual(values = fill_genomicCompartment)

#pdf("DegenerationTempoWith_meanFOP_vs_meandS_withSE.pdf")
#	p1_withSE
#dev.off()


#pdf("DegenerationTempoWith_meanFOP_vs_meandS_withoutSE.pdf")
#	p1_withoutSE
#dev.off()



