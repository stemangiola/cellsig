# Process data

library(tidyverse)
library(cellsig)
library(tidybulk)

# 3 Dataset GSE135390_raw_counts.csv
GSE135390_Th22 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th22=lib413, BC120926_Th22=lib424, BC121206_Th22=lib435) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper", level=4,
         note = "CD127+CD25-CD45RA-CXCR3-CCR6+CCR4+CCR10+ Th22 sorted from peripheral blood")


GSE135390_Th1_17 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th1_17 = lib415, BC120926_Th1_17 = lib426, BC121206_Th1_17 = lib437) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper", level=4,
         note = "CD127+CD25-CD45RA-CXCR3+CCR6+ Th1_17 sorted from peripheral blood")

GSE135390_new <- bind_rows(GSE135390_Th22, GSE135390_Th1_17) %>% mutate(count=as.numeric(count))


# 11 Dataset GSE115103_raw_counts
######### NEW 
GSE115103_Th1_17 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165172=S3858, GSM3165182=S3871, GSM3165192=S3885) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper", level=4,
         note="CD4+ Th1_Th17 cells sorted from PBMC")

GSE115103_Th22 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165173=S3859, GSM3165183=S3872, GSM3165193=S3886) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper", level=4,
         note="CD4+ Th22 cells sorted from PBMC")

GSE115103_new <- bind_rows(GSE115103_Th1_17, GSE115103_Th22)

# 20 Dataset GSE123812_CD4_bulk_RNA_counts
########## New

GSE123812_TFH1_17 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, GSM3511701=`Tfh1-17_1`, GSM3511702=`Tfh1-17_2`, GSM3511703=`Tfh1-17_3`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper", level=4,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3+CCR6+ CD4+ primary follicular Th1-17 cells sorted from PBMCs")

GSE123812_TFH1 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511704=Tfh1_1, GSM3511705=Tfh1_2, GSM3511706=Tfh1_3, GSM3511707=Tfh1_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper", level=4,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3+CCR6- CD4+ follicular Th1 cells sorted from PBMCs")

GSE123812_TFH17 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511708=Tfh17_1, GSM3511709=Tfh17_2, GSM3511710=Tfh17_3, GSM3511711=Tfh17_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper", level=4,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3-CCR6+ CD4+ follicular Th17 cells sorted from PBMCs")


GSE123812_TFH2 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511712=Tfh2_1, GSM3511713=Tfh2_2, GSM3511714=Tfh2_3, GSM3511715=Tfh2_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper", level=4,
         note="CXCR5+ CD4+CD25-IL7RhiCD45RA-CXCR3-CCR6- CD4+ follicular Th2 cells sorted from PBMCs")

GSE123812_TH1_17 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511716=`Th1-17_1`, GSM3511717=`Th1-17_2`, GSM3511718=`Th1-17_3`, GSM3511719=`Th1-17_4`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper", level=4,
         note="CD25-IL7RhiCD45RA-CXCR3+CCR6+ CD4+ Th1-17 cells sorted from PBMCs")

GSE123812_new <- bind_rows(GSE123812_TFH1_17, GSE123812_TFH1, GSE123812_TFH17, 
                           GSE123812_TFH2, GSE123812_TH1_17)


# below are data collected for new_data3
# 26  GSE164643 Dataset GSE164643_Raw_counts_moDC_all_samples_txt
GSE164643_moDC_immature <- GSE164643_Raw_counts_moDC_all_samples_txt %>% 
  select(symbol = GENE_ID, GSM5016287 = A1, GSM5016288 = D1stMED, GSM5016289 = D2, 
         GSM5016290 = D2stMED, GSM5016291	 = D3stMED, GSM5016292 = D6stMED) %>% 
  pivot_longer(-symbol, names_to = "sample", values_to = "count") %>% 
  mutate(dataset = "GSE164643", cell_type = "dendritic_myeloid_immature", level=4,
         note="monocyte-derived DCs derived from blood CD14+ monocytes cultured for 6d with GM-CSF IL-4")

GSE164643_moDC_mature <- GSE164643_Raw_counts_moDC_all_samples_txt %>% 
  select(symbol = GENE_ID, GSM5016233 = B1, GSM5016234=D1stA, GSM5016235=D2stA, GSM5016236=D3stA,
         GSM5016237=D6stA, GSM5016238=E2, GSM5016239=C1, GSM5016240=D1stB, GSM5016241=D2stB,
         GSM5016242=D3stB, GSM5016243=D6stB, GSM5016244=F2, GSM5016245=D1, GSM5016246=D1stD,
         GSM5016247=D2stD, GSM5016248=D3stD, GSM5016249=D6stD, GSM5016250=G2, GSM5016251=D1stF,
         GSM5016252=D2stF, GSM5016253=D3stF, GSM5016254=D6stF, GSM5016255=E1, GSM5016256=H2,
         GSM5016257=B3, GSM5016258=D1stG, GSM5016259=D2stG, GSM5016260=D3stG, GSM5016261=D6stG,
         GSM5016262=F1, GSM5016263=C3, GSM5016264=D1stH, GSM5016265=D2stH, GSM5016266=D3stH,
         GSM5016267=D6stH, GSM5016268=G1, GSM5016269=D1stJ, GSM5016270=D2stJ, GSM5016271=D3,
         GSM5016272=D3stJ, GSM5016273=D6stJ, GSM5016274=H1, GSM5016275=A2, GSM5016276=A3, 
         GSM5016277=D1stK, GSM5016278=D2stK, GSM5016279=D3stK, GSM5016280=D6stK, GSM5016281=B2,
         GSM5016282=D1stL, GSM5016283=D2stL, GSM5016284=D3stL, GSM5016285=D4stL, GSM5016286=D6stL) %>% 
  pivot_longer(-symbol, names_to = "sample", values_to = "count") %>% 
  mutate(dataset = "GSE164643", cell_type = "dendritic_myeloid_mature", level=4,
         note="monocyte-derived DCs activated by various strains of B.pertussis for 3hrs")

GSE164643 <- bind_rows(GSE164643_moDC_immature, GSE164643_moDC_mature)

# 27 GSE60424 TAb delimited file containing TMM normalised counts
GSE60424_neutrophil <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib226, GSM1479499=lib288, GSM1479506=lib295, GSM1479520=lib309) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="neutrophil", level=3,
         note="neutrophil from PBMC")

GSE60424_monocyte <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib227, GSM1479500=lib289, GSM1479507=lib296, GSM1479521=lib310) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="monocyte", level=3,
         note="monocyte from PBMC")

GSE60424_B_cell <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib228, GSM1479501=lib290, GSM1479508=lib297, GSM1479522=lib311) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="b_cell", level=2,
         note="b cell from PBMC")

GSE60424_CD4 <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib229, GSM1479502=lib291, GSM1479509=lib298, GSM1479523=lib312) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="t_CD4", level=3,
         note="t_CD4 cell from PBMC")

GSE60424_CD8 <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib230, GSM1479503=lib292, GSM1479510=lib299, GSM1479524=lib313) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="t_CD8", level=3,
         note="t_CD8 cell from PBMC")

GSE60424_NK <- GSE60424_GEOSubmit_FC1to11_normalized_counts_txt %>% 
  select(ensembl=genenames, GSM1479438=lib231, GSM1479504=lib293, GSM1479511=lib300, GSM1479525=lib314) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  rename(symbol=transcript) %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE60424", cell_type="nk_resting", level=3,
         note="NK cell from PBMC of healthy individual")

GSE60424 <- bind_rows(GSE60424_neutrophil, GSE60424_monocyte, GSE60424_B_cell,
                      GSE60424_CD4, GSE60424_CD8, GSE60424_NK)

# 28 GSE107011 Gianni Monarco

geneNames <- kallisto_gene$genes %>% 
  as_tibble() %>% 
  select(gene_id, gene_name)

GSE107011_data <- kallisto_gene$counts %>% 
  as_tibble(rownames = "ensembl") %>% 
  mutate(across(where(is.double), ceiling)) %>% 
  left_join(geneNames, by = c("ensembl"="gene_id")) %>% 
  rename(symbol = gene_name) %>% 
  mutate(ensembl = str_extract(ensembl, ".*(?=\\.)"))
# select(-ref_genome) %>% 
# rename(symbol=transcript)

GSE107011_CD8_naive <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859411=RHH5201, GSM2859439=RHH5229, GSM2859468=RHH5258, GSM2859497=RHH5287) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD8", level=3,
         note="t_CD8_naive cell from PBMC of healthy individual")

GSE107011_CD8_CM <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859412=RHH5202, GSM2859440=RHH5230, GSM2859469=RHH5259, GSM2859498=RHH5288) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD8_memory_central", level=5,
         note="t_CD8 central memory cell from PBMC of healthy individual")

GSE107011_CD8_EM <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859413=RHH5203, GSM2859441=RHH5231, GSM2859470=RHH5260, GSM2859499=RHH5289) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD8_memory_effector", level=5,
         note="t_CD8 effector memory cell from PBMC of healthy individual")

GSE107011_CD8_TE <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859414=RHH5204, GSM2859442=RHH5232, GSM2859471=RHH5261, GSM2859506=RHH5296) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD8", level=3,
         note="t_CD8 terminal effector cell from PBMC of healthy individual")

GSE107011_MAIT <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859415=RHH5205, GSM2859443=RHH5233, GSM2859472=RHH5262, GSM2859507=RHH5297) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_cell", level=2,
         note="MAIT cell from PBMC of healthy individual")

GSE107011_gdT_VD2pos <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859416=RHH5206, GSM2859444=RHH5234, GSM2859473=RHH5263, GSM2859508=RHH5298) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_gamma_delta", level=3,
         note="VD2+ gamma delta t cell from PBMC of healthy individual")

GSE107011_gdT_VD2neg <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859417=	RHH5207,  GSM2859445=RHH5235, GSM2859474=RHH5264, GSM2859509=RHH5299) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_gamma_delta", level=3,
         note="VD2- gamma delta t cell from PBMC of healthy individual")

GSE107011_TFH <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859418=RHH5208, GSM2859446=RHH5236, GSM2859475=RHH5265, GSM2859510=RHH5300) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_helper", level=4,
         note="follicular helper t cell from PBMC of healthy individual")

GSE107011_Treg <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859419=RHH5209, GSM2859447=RHH5237, GSM2859476=RHH5266, GSM2859511=RHH5301) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_reg", level=5,
         note="regulatory t cell from PBMC of healthy individual")

GSE107011_Th1 <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859420=RHH5210, GSM2859448=RHH5238, GSM2859477=RHH5267, GSM2859512=RHH5302) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_helper_h1", level=5,
         note="Th1 cell from PBMC of healthy individual")

GSE107011_Th1_Th17 <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859421=RHH5211, GSM2859449=RHH5239, GSM2859478=RHH5268, GSM2859513=RHH5303) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_helper", level=4,
         note="Th1/Th17 cell from PBMC of healthy individual")

GSE107011_Th17 <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859422=RHH5212, GSM2859450=RHH5240, GSM2859479=RHH5269, GSM2859514=RHH5304) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_helper_h17", level=5,
         note="Th17 cell from PBMC of healthy individual")

GSE107011_Th2 <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859423=RHH5213, GSM2859451=RHH5241, GSM2859480=RHH5270, GSM2859515=RHH5305) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_helper_h2", level=5,
         note="Th2 cell from PBMC of healthy individual")

GSE107011_t_CD4_naive <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859424=RHH5214, GSM2859452=RHH5242, GSM2859481=RHH5271, GSM2859516=RHH5306) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD4", level=3,
         note="Naive CD4 t cell from PBMC of healthy individual")

GSE107011_t_CD4_TE <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859453=RHH5243, GSM2859482=RHH5272) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="t_CD4", level=3,
         note="Terminal Effector CD4 t cell from PBMC of healthy individual")

# do not include this progenitor cell data
# GSE107011_progenitor <- GSE107011_data %>% 
#   select(ensembl, symbol, GSM2859425=RHH5215, GSM2859454=RHH5244, GSM2859483=RHH5273, GSM2859517=RHH5307) %>% 
#   pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
#   mutate(dataset="GSE107011", cell_type="progenitor", level=?,
#          note="progenitor cell from PBMC of healthy individual")

GSE107011_b_naive <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859426=RHH5216, GSM2859455=RHH5245, GSM2859484=RHH5274, GSM2859518=RHH5308) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="b_naive", level=3,
         note="Naive b cell from PBMC of healthy individual")

GSE107011_b_NSM <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859427=RHH5217, GSM2859456=RHH5246, GSM2859485=RHH5275, GSM2859519=RHH5309) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="b_memory", level=3,
         note="Non-switched memory b cell from PBMC of healthy individual")

GSE107011_ExB <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859428=RHH5218, GSM2859457=RHH5247, GSM2859486=RHH5276, GSM2859520=RHH5310) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="b_memory", level=3,
         note="Exhausted b cell from PBMC of healthy individual")

GSE107011_B_SM <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859429=RHH5219, GSM2859458=RHH5248, GSM2859487=RHH5277, GSM2859521=RHH5311) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="b_memory", level=3,
         note="Switched memory b cell from PBMC of healthy individual")

GSE107011_plasmablast <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859430=RHH5220, GSM2859459=RHH5249, GSM2859488=RHH5278, GSM2859522=RHH5312) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="b_cell", level=2,
         note="Plasmablasts from PBMC of healthy individual")

GSE107011_classic_Mono <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859431=RHH5221, GSM2859460=RHH5250, GSM2859489=RHH5279, GSM2859523=RHH5313) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="monocyte", level=3,
         note="Classical monocytes from PBMC of healthy individual")

GSE107011_intermediate_Mono <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859432=RHH5222, GSM2859461=RHH5251, GSM2859490=RHH5280, GSM2859524=RHH5314) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="monocyte", level=3,
         note="Intermediate monocytes from PBMC of healthy individual")

GSE107011_nonClassic_Mono <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859433=RHH5223, GSM2859462=RHH5252, GSM2859491=RHH5281, GSM2859525=RHH5315) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="monocyte", level=3,
         note="Non-classical monocytes from PBMC of healthy individual")

GSE107011_NK <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859434=RHH5224, GSM2859463=RHH5253, GSM2859492=RHH5282, GSM2859526=RHH5316) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="nk_resting", level=3,
         note="Natural killer cell from PBMC of healthy individual")

# do not include this plasmacytoid dendritic cell data
# GSE107011_pDC <- GSE107011_data %>% 
#   select(ensembl, symbol, GSM2859435=RHH5225, GSM2859464=RHH5254, GSM2859493=RHH5283, GSM2859527=RHH5317) %>% 
#   pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
#   mutate(dataset="GSE107011", cell_type="dentritic_plsamacytoid", level=?,
#          note="Plasmacytoid DC from PBMC of healthy individual")

GSE107011_mDC <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859436=RHH5226, GSM2859465=RHH5255, GSM2859494=RHH5284, GSM2859528=RHH5318) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="dentritic_myeloid", level=3,
         note="Myeloid DC from PBMC of healthy individual")

GSE107011_neutrophil <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859437=RHH5227, GSM2859466=RHH5256, GSM2859495=RHH5285, GSM2859529=RHH5319) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="neutrophil", level=3,
         note="Neutrophil from PBMC of healthy individual")

GSE107011_basophil <- GSE107011_data %>% 
  select(ensembl, symbol, GSM2859438=RHH5228, GSM2859467=RHH5257, GSM2859496=RHH5286, GSM2859530=RHH5320) %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE107011", cell_type="granulocyte", level=2,
         note="Basophil from PBMC of healthy individual")

GSE107011 <- bind_rows(
  GSE107011_CD8_naive, GSE107011_CD8_CM, GSE107011_CD8_EM, GSE107011_CD8_TE, GSE107011_MAIT, 
  GSE107011_gdT_VD2pos, GSE107011_gdT_VD2neg, GSE107011_TFH, GSE107011_Treg, GSE107011_Th1, GSE107011_Th1_Th17, 
  GSE107011_Th17, GSE107011_Th2, GSE107011_t_CD4_naive, GSE107011_t_CD4_TE, GSE107011_b_naive, GSE107011_b_NSM,
  GSE107011_ExB, GSE107011_B_SM, GSE107011_plasmablast, GSE107011_classic_Mono, GSE107011_intermediate_Mono, 
  GSE107011_nonClassic_Mono, GSE107011_NK, GSE107011_mDC, GSE107011_neutrophil, GSE107011_basophil
)

# 29 GSE157844

GSE157844_DC_immature1 <- GSM4776603_I20_1161_01_htseq_txt %>% 
  select(ensembl = X1, symbol=X2, count = X3) %>% 
  mutate(dataset="GSE157844", sample="GSM4776603", cell_type = "dendritic_myeloid_immature", level=4,
         note="monocyte derived DCs: monocyte from PBMC cultured 7d with GM-CSF IL-4")

GSE157844_DC_immature2 <- GSM4776604_I20_1161_03_htseq_txt %>% 
  select(ensembl = X1, symbol=X2, count = X3) %>% 
  mutate(dataset="GSE157844", sample="GSM4776604", cell_type = "dendritic_myeloid_immature", level=4,
         note="monocyte derived DCs: monocyte from PBMC cultured 7d with GM-CSF IL-4")

GSE157844_DC_mature1 <- GSM4776601_I20_1161_02_htseq_txt %>% 
  select(ensembl = X1, symbol=X2, count = X3) %>% 
  mutate(dataset="GSE157844", sample="GSM4776601", cell_type = "dendritic_myeloid_mature", level=4,
         note="monocyte derived DCs stimulated for 3d by LPS")

GSE157844_DC_mature2 <- GSM4776601_I20_1161_02_htseq_txt %>% 
  select(ensembl = X1, symbol=X2, count = X3) %>% 
  mutate(dataset="GSE157844", sample="GSM4776602", cell_type = "dendritic_myeloid_mature", level=4,
         note="monocyte derived DCs stimulated for 3d by LPS")

GSE157844 <- bind_rows(GSE157844_DC_immature1, GSE157844_DC_immature2,
                       GSE157844_DC_mature1, GSE157844_DC_mature2)

# 30 GSE174659 CD141+ DC 

# CD141+ cDC1
GSE174659_cDC1 <- GSE174659_Raw_counts_csv %>% 
  select(symbol=X1, GSM5321871=HB141_02, GSM5321872=HB141_16, GSM5321873=HB141_24, GSM5321874=HB141_25,
         GSM5321875=HB141_27, GSM5321876=HB141_28, GSM5321877=HB141_29, GSM5321913=HP141_02, 
         GSM5321914=HP141_16, GSM5321915=HP141_24, GSM5321916=HP141_25, GSM5321917=HP141_26, 
         GSM5321918=HP141_27, GSM5321919=HP141_28, GSM5321920=HP141_29) %>% 
  pivot_longer(-symbol, names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE174659", cell_type = "dendritic_myeloid", level=3, 
         note="CD141+ cDC1 isolated from bronchoalveolar lavage of healthy controls")

# CD1c+ cDC2
GSE174659_cDC2 <- GSE174659_Raw_counts_csv %>% 
  select(symbol=X1, GSM5321878=HB1C_02, GSM5321879=	HB1C_16, GSM5321880=HB1C_24, GSM5321881=HB1C_25,
         GSM5321882=HB1C_26, GSM5321883=HB1C_27, GSM5321884=HB1C_28, GSM5321885=HB1C_29, GSM5321886=HB1C_55,
         GSM5321921=HP1C_02, GSM5321922=HP1C_16, GSM5321923=HP1C_24, GSM5321924=HP1C_25, GSM5321925=HP1C_26,
         GSM5321926=HP1C_27, GSM5321927=HP1C_28, GSM5321928=HP1C_29, GSM5321929=HP1C_55) %>% 
  pivot_longer(-symbol, names_to="sample", values_to="count") %>% 
  mutate(dataset="GSE174659", cell_type = "dendritic_myeloid", level=3, 
         note="CD1+ cDC2 isolated from bronchoalveolar lavage of healthy controls")

GSE174659 <- bind_rows(GSE174659_cDC1, GSE174659_cDC2)

# 31 

counts_third_db_raw <- bind_rows(GSE164643, GSE60424, GSE157844, GSE174659, GSE107011, 
                       GSE135390_new, GSE115103_new, GSE123812_new)


saveRDS(counts_third_db_raw, file = "dev/counts_third_db_raw.rds")


