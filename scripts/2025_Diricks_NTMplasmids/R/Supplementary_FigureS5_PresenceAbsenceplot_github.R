#Load libraries
library(tidyverse)
library(cowplot)
library(janitor)
library(readxl)

#Read in data-------------------------------------------------------------------------------------------------
#Read in fastANI results
fast_ani_multiple <- read_excel("fastANI_mycobacteriaceaeAssemblies_results_v2.xlsx",
                                sheet="Summary >1 species")
colnames(fast_ani_multiple) <- make_clean_names(colnames(fast_ani_multiple))
View(fast_ani_multiple)

# Create graph for plasmids detected in more than one species-------------------------------
#Reorder species by branch order from phylogenetic tree and introduce missing species
fast_ani_multiple$species <- factor(fast_ani_multiple$species, levels  = c("parascrofulaceum", "scrofulaceum", "nebraskense", "europaeum", "seoulense", "paraseoulense", "paraffinicum", "helveticum", "parmense", "vulneris", "senriense", "colombiense", "mantenii", "intracellulare", "paraintracellulare", "marseillense", "avium", "interjectum", "terramassiliense", "paraense", "alsense", "heidelbergense", "malmoense", "palustre", "conspicuum", "bohemicum", "cambodiense", "saskatchewanense", "florentinum", "stomatepiae", "triplex", "lentiflavum", "montefiorense", "simiae", "sherrisii", "numidiamassiliense", "rhizamassiliense", "shigaense", "persicum", "innocens", "kansasii", "pseudokansasii", "gastri", "ostraviense", "attenuatum", "ulcerans","liflandii", "pseudoshottsii", "marinum", "shottsii", "basiliense", "simulans", "riyadhense", "lacus", "decipiens", "spongiae", "shinjukuense", "haemophilum", "szulgai", "angelicum", "bourgelatii", "intermedium", "asiaticum", "gordonae", "paragordonae", "vicinigordonae", "kiyosense", "kubicae", "celatum", "branderi", "kyorinense", "fragae", "shimoidei", "heckeshornense", "xenopi", "noviomagense", "botniense", "cookii", "paraterrae", "senuensis", "algericum", "kumamotonensis", "terrae", "longobardus", "nativiensis", "crassicus", "vasticus", "zoologicum", "virginiensis", "icosiumassiliensis", "heraklionensis", "nonchromogenicus", "arupensis", "minnesotensis", "hiberniae", "acidiphilus", "koreensis", "parakoreensis", "trivialis", "talmoniae", "fallax", "insubricum", "brumae", "palauense", "rhodesiae", "aichiense", "aromaticivorans", "pallens", "crocinum", "sarraceniae", "helvum", "sphagni", "vinylchloridicum", "pinniadriaticum", "anyangense", "alkanivorans", "malmesburyense", "komanii", "novocastrense", "neumannii", "celeriflavum", "rutilum", "elephantis", "pulveris", "holsaticum", "deserti", "tusciae", "gadium", "neglectum", "moriokaense", "barrassiae", "stellerae", "agri", "hubeiense", "thermoresistibile", "hassiacum", "phlei", "doricum", "monacense", "litorale", "baixiangningiae", "gossypii", "manitobense", "vanbaalenii", "vaccae", "parafortuitum", "iranicum", "bulgaricum", "aurum", "hippocampi", "poriferae", "gilvum","austroafricanum", "pyrenivorans", "duvalii", "psychrotolerans", "chubuense","chlorophenolicum", "obuense", "xanthum", "senegalense","conceptionense","farcinogenes", "syngnathidarum", "boenickei", "porcinum", "peregrinum", "fortuitum", "setense", "alvei", "lutetiense", "houstonense", "fortunisiensis", "aquaticum", "dioxanotrophicus", "brisbanense", "mageritense", "wolinskyi", "smegmatis", "goodii", "arabiense", "lacusdiani", "sediminis", "grossiae", "arenosum", "hodleri", "yunnanensis", "madagascariense", "chelonae", "stephanolepidis", "salmoniphilum", "saopaulense", "immunogenum", "franklinii", "abscessus", "phocaicum", "llatzerense", "aubagnense", "NA", "neoaurum.II", "neoaurum.I", "bacteremicum", "adipatum", "diernhoferi", "cosmeticum", "canariasense", "fluranthenivorans", "tokaiense", "mengxianglii", "komossense", "confluentis", "chitae","Novel species"))

#Filter out novel species
fast_ani_multiple <- filter(fast_ani_multiple,species!="Novel species")

#Bubble plot with plasmids
dot_plot <- ggplot(fast_ani_multiple)+
  geom_tile(aes(y=species,
                x=reorder(plasmid,plasmid,function(x)-length(x)),
                fill=growth))+
  theme(#axis.text.x = element_text(angle=90, hjust=1)
  axis.text.x=element_blank()
    )+
  labs(x="",y="Species in which plasmid was detected",title="A",fill="Runyon-classification")+
  scale_fill_manual(values=c("#d2691e","#92e2fa","#1602f7","black","darkgrey"))
dot_plot

#Bar plot for plasmid detection by growth type - Number of strains
bar_plot <- ggplot(fast_ani_multiple)+
  geom_bar(aes(x = reorder(plasmid,plasmid,function(x)-length(x)),
               fill=growth))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(x="Plasmid",y="Number of strains",fill="Runyon-classification",title="B")+
  scale_fill_manual(values=c("#d2691e","#92e2fa","#1602f7","black","darkgrey"))
bar_plot

#Put graphs together
plot_grid(dot_plot,bar_plot,ncol=1,align = c("hv"),
          axis = c("tblr"), rel_heights = c(2,0.8))