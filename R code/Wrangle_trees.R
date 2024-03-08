#Bethany Allen (27th May 2022)
#Code to extract tip labels for defining clades

library(tidyverse)
library(ape)
library(palaeoverse)


#Sample trees from output of Lloyd et al. (2016) to create MCC tree

trees <- read.tree("trees/Lloyd2_dinosaur_str_1000_mpts.tre")
sample_numbers <- sample(x = length(trees), size = 1000, replace = F)
tree_subsample <- trees[c(sample_numbers)]
write.tree(tree_subsample, "trees/Lloyd2_sample1000.tree")


#Check trees against PBDB data

tree1 <- read.tree("trees/Lloyd1.tree")
tree2 <- read.tree("trees/Benson1.tree")
tree3 <- read.tree("trees/Lloyd2_MCC.tree")
no_NA <- read_csv("trees/tip_constraints.csv")

test <- phylo_check(tree = tree1, list = no_NA$taxon_name, out = "diff_table",
                    sort = "presence")
test <- filter(test, present_in_tree == TRUE)


#Remove birds, informal taxa and invalid taxa
#See 'Tree modifications log' for list of changes made

tree1 <- read.tree("trees/Lloyd1.tree")

Lloyd1_to_remove <- c("Chialingosaurus_kuani", "Lexovisaurus_durobrivensis",
                      "Stegosaurus_armatus", "Minmi_paravertebra",
                      "Tianzhenosaurus_youngi", "Stormbergia_dangershoeki",
                      "Ornatotholus_browni", "Dwarf_pachycephalosaur",
                      "Stygimoloch_spinifer", "Psittacosaurus_guyangensis",
                      "Nedoceratops_hatcheri", "Fontllonga_euhadrosaur",
                      "Ammosaurus_major", "Apatosaurus_sp.",
                      "Pleurocoelus_nanus", "Titanosaurus_indicus",
                      "Japalpur_titanosaur_indet.", "Malagasy_Taxon_B",
                      "Toba_titanosaur", "Peiropolis_titanosaur",
                      "Bajo_Barreal_abelisaur", "Lameta_abelisaur",
                      "Cristatusaurus_lapparenti",
                      "Two_Medicine_tyrannosaurine", "Quarry_nine_coelurosaur",
                      "El_Brete_oviraptor", "Kamyn_Khondt_oviraptorine",
                      "Gugyedong_maniraptoran", "Dornogov_troodontid")

tree1 <- drop.tip(tree1, Lloyd1_to_remove)
write.tree(tree1, "trees/Lloyd1clean.tree")

tree2 <- read.tree("trees/Benson1.tree")

Benson_to_remove <- c("Minmi_paravertebra", "Zhongyuanosaurus_luoyangensis",
                      "Crichtonsaurus_bohlini", "Minotaurasaurus_ramachandrani",
                      "Tianzhenosaurus_youngi", "Antarctopelta_oliveroi",
                      "Zhejiangosaurus_lishuiensis", "Stormbergia_dangershoeki",
                      "Nedoceratops_hatcheri", "Willinakaqe_salitralensis",
                      "Kundurosaurus_nagornyi", "Sahaliyania_elunchunorum",
                      "Chuxiongosaurus", "French_Bothriospondylus",
                      "Aeolosaurus_sp.", "CV00214", "IGM100/1323",
                      "IGM100/1126", "IGM100/44", "Sapeornis_chaoyangensis",
                      "Jeholornis_prima", "Jixiangornis_orientalis",
                      "Eoconfuciusornis_zhengi",
                      "Changchengornis_hengdaoziensis", "Confuciusornis_sanctus",
                  "Confuciusornis_dui", "Jinzhouornis_zhangjiyingia",
                  "Protopteryx_fengningensis", "Otogornis_genghisi",
                  "Elsornis_keni", "Shenqiornis_mengi",
                  "Longipteryx_chaoyangensis", "Boluochia_zhengi",
                  "Rapaxavis_pani", "Iberomesornis_romerali",
                  "Shanweiniao_cooperorum", "Longirostravis_hani",
                  "Vescornis_hebeiensis", "Pengornis_houi",
                  "Gobipteryx_minuta", "Neuquenornis_volans",
                  "Eoenantiornis_buhleri", "Concornis_lacustris",
                  "Eocathayornis_walkeri", "Cathayornis_yandica",
                  "Liaoningornis_longidigitris", "Eoalulavis_hoyasi",
                  "Archaeorhynchus_spathula", "Patagopteryx_deferrariisi",
                  "Jianchangornis_microdonta", "Schizooura_lii",
                  "Vorona_berivotrensis", "Zhongjianornis_yangi",
                  "Chaoyangia_beishanensis", "Hongshanornis_longicresta",
                  "Longicrusavis_houi", "Yixianornis_grabaui",
                  "Yanornis_martini", "Songlingornis_linghensis",
                  "Gansus_yumenensis", "Apsaravis_ukhaana",
                  "Ambiortus_dementjevi", "Hollanda_luceria",
                  "Ichthyornis_dispar", "Vegavis_iaai", "Limenavis_patagonica",
                  "Enaliornis_barretti", "Baptornis_advenus",
                  "Baptornis_varneri", "Parahesperornis_alexi",
                  "Hesperornis_regalis", "Dromiceomimus_brevitertius")

tree2 <- drop.tip(tree2, Benson_to_remove)
write.tree(tree2, "trees/Benson1clean.tree")

tree3 <- read.tree("trees/Lloyd2_MCC.tree")

Lloyd2_to_remove <- c(grep("^\\D+$", tree3$tip.label, value = TRUE, invert = TRUE),
                      "allzero", "Sapeornis_chaoyangensis",
                      "Jeholornis_prima", "Jixiangornis_orientalis",
                      "Eoconfuciusornis_zhengi",
                      "Changchengornis_hengdaoziensis", "Confuciusornis_sanctus",
                      "Confuciusornis_dui", "Jinzhouornis_zhangjiyingia",
                      "Protopteryx_fengningensis", "Otogornis_genghisi",
                      "Elsornis_keni", "Shenqiornis_mengi",
                      "Longipteryx_chaoyangensis", "Boluochia_zhengi",
                      "Rapaxavis_pani", "Iberomesornis_romerali",
                      "Shanweiniao_cooperorum", "Longirostravis_hani",
                      "Vescornis_hebeiensis", "Pengornis_houi",
                      "Gobipteryx_minuta", "Neuquenornis_volans",
                      "Eoenantiornis_buhleri", "Concornis_lacustris",
                      "Eocathayornis_walkeri", "Cathayornis_yandica",
                      "Liaoningornis_longidigitris", "Eoalulavis_hoyasi",
                      "Archaeorhynchus_spathula", "Patagopteryx_deferrariisi",
                      "Jianchangornis_microdonta", "Schizooura_lii",
                      "Vorona_berivotrensis", "Zhongjianornis_yangi",
                      "Chaoyangia_beishanensis", "Hongshanornis_longicresta",
                      "Longicrusavis_houi", "Yixianornis_grabaui",
                      "Yanornis_martini", "Songlingornis_linghensis",
                      "Gansus_yumenensis", "Apsaravis_ukhaana",
                      "Ambiortus_dementjevi", "Hollanda_luceria",
                      "Ichthyornis_dispar", "Vegavis_iaai", "Limenavis_patagonica",
                      "Enaliornis_barretti", "Baptornis_advenus",
                      "Baptornis_varneri", "Parahesperornis_alexi",
                      "Hesperornis_regalis", "Dromiceomimus_brevitertius",
                      "Cimolopteryx_minima", "Cimolopteryx_rara",
                      "Cimolopteryx_petra", "Austinornis_lentus",
                      "Ceramornis_major", "Cimolopteryx_maxima",
                      "Palintropus_retusus", "Polarornis_gregorii",
                      "Torotix_clemensi", "Telmatornis_priscus",
                      "Apatornis_celer", "Guildavis_tener", "Iaceornis_marshi",
                      "Brodavis_americanus", "Yanornis_wui",
                      "Avisaurus_archibaldi", "Avisaurus_gloriae",
                      "Soroavisaurus_australis", "Huoshanornis_huji",
                      "Mystiornis_cyrili", "Halimornis_thompsoni",
                      "Enantiophoenix_electrophyla", "Enantiornis_leali",
                      "Sinornis_santensis", "Xiangornis_shenmi",
                      "Noguerornis_gonzalezi", "Nanantius_eos",
                      "Yungavolucris_brevipedalis", "Qiliania_graffini",
                      "Sulcavis_geeorum", "Bohaiornis_guoi",
                      "Alexornis_antecedens", "Dalingheornis_liweii",
                      "Lectavis_bretincola", "Shenzhouraptor_sinensis",
                      "Brachiosauridae_indet_MHNM_Unreferred",
                      "Nomingia_gobiensis", "Sigilmassasaurus_brevicollis",
                      "Sauroniops_pachytholus", "Trigonosaurus_pricei",
                      "Paluxysaurus_jonesi", "Brontomerus_mcintoshi",
                      "Astrodon_johnstoni", "Chuxiongosaurus_lufengensis",
                      "Kundurosaurus_nagornyi", "Willinakaqe_salitralensis",
                      "Sahaliyania_elunchunorum", "Kukufeldia_tilgatensis",
                      "Mojoceratops_perifania", "Stygimoloch_spinifer",
                      "Texacephale_langstoni", "Stormbergia_dangershoeki",
                      "Tianzhenosaurus_youngi", "Zhongyuanosaurus_luoyangensis",
                      "Minmi_paravertebrata", "Antarctopelta_oliveroi",
                      "Zhejiangosaurus_lishuiensis", "Stegosaurus_armatus",
                      "Uteodon_aphanoecetes", "Zhongyuangosaurus_luoyangensis",
                      "Dollodon_bampingi", "Prosaurolophus_blackfeetensis",
                      "Anatotitan_copei", "Edmontosaurus_saskatchewanensis")

tree3 <- drop.tip(tree3, Lloyd2_to_remove)
write.tree(tree3, "trees/Lloyd2clean.tree")


#Create subtrees

#Read in tree
tree <- read.tree("trees/Lloyd1clean.tree")

#Extract tip labels
species_list <- tree$tip.label

#Write text file of tip labels
species_file <- file("trees/Lloyd1_clades.txt")
writeLines(species_list, species_file)
close(species_file)

#Manually split tip names into major clades

#Read in subclade list
clade <- read.table("trees/Lloyd1_Theropoda.txt")
clade <- clade$V1

#Check that subclade clusters on phylogeny
plot(tree, show.tip.label = F)
tiplabels(tip = match(clade, tree$tip.label), pch = 19, col = "red")

#Trim phylogeny to tips in subclade
clade_tree <- phylo_check(tree = tree, list = clade, out = "tree")

#Save subclade phylogeny
write.tree(clade_tree, "trees/Lloyd1_Theropoda.tree")


#Plot the distribution of age ranges over time
no_NA <- read_csv("trees/tip_constraints.csv")
tree <- read.tree("trees/Lloyd2clean.tree")

tips_in_trees <- filter(no_NA, taxon_name %in% tree$tip.label)
tips_in_trees <- tips_in_trees[order(tips_in_trees$lastapp_min_ma), ]

ggplot(data = tips_in_trees, aes(xmin = lastapp_min_ma,
                                 xmax = firstapp_max_ma,
                                 y = seq(1, nrow(tips_in_trees)),
                                 group = 1)) +
  annotate(geom = "rect", xmin = c(23.8, 79.0, 108.7, 171.0),
           xmax = c(0, 55.4, 95.5, 135.4),
           ymin = -Inf, ymax = Inf,
           colour = "grey", alpha = 0.2, linewidth = 0)  +
  geom_errorbar(linewidth = 1, width = 0.8) +
  scale_x_reverse() +
  labs(x = "Time (My before K-Pg boundary)", y = "Taxon ID") +
  theme_classic(base_size = 14)
ggsave("figures/Sampling_through_time.pdf", scale = 3.5, dpi = 600)
