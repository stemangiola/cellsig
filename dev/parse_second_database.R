library(tidybulk)
library(tidyverse)
library(ggplot2)


# 1 Dataset GSE151586_RAW # DECIMAL COUNTS!
GSM4586306_aTreg <- GSM4586306_S143d3IL2_abundance_tsv  %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586306", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586309_aTreg <- GSM4586309_S146d3IL2_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586309", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586312_aTreg <- GSM4586312_HA5231IL2_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586312", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586315_aTreg <- GSM4586315_HA5232IL2_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586315", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586318_aTreg <- GSM4586318_HA5233IL2_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586318", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586321_aTreg <- GSM4586321_HA5526_C_abundance_tsv %>%  
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586321", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586324_aTreg <- GSM4586324_HA5527_C_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586324", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586327_aTreg <- GSM4586327_HA5528_C_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586327", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSM4586330_aTreg <- GSM4586330_HA5529_C_abundance_tsv %>% 
  select(ensembl = target_id, count = est_counts) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  mutate(sample = "GSM4586330", dataset = "GSE151586", cell_type="t_reg", level=5,
         note = "CD4+CD25hiCD127lo/neg Tregs sorted from peripheral blood IL2/antiCD3/anti_CD28_activated for 3 days")

GSE151586 <- bind_rows(GSM4586306_aTreg, GSM4586309_aTreg, GSM4586312_aTreg, GSM4586315_aTreg,
                       GSM4586318_aTreg, GSM4586321_aTreg, GSM4586324_aTreg, GSM4586327_aTreg,
                       GSM4586330_aTreg)


# 2 Dataset GSE151536_RAW
GSM4584330_Treg <- GSM4584330_Blank1_txt %>% select(ensembl = id, symbol = Symbol, count = Blank1_count) %>% 
  mutate(sample = "GSM4584330", dataset = "GSE151536", cell_type="t_reg", level=5,
         note = "Tregs sorted from PBMC") 

GSM4584332_Treg <- GSM4584332_Blank2_txt %>% select(ensembl = id, symbol = Symbol, count = Blank2_count) %>% 
  mutate(sample = "GSM4584332", dataset = "GSE151536", cell_type="t_reg", level=5,
         note = "Tregs sorted from PBMC")

GSM4584334_Treg <- GSM4584334_Blank3_txt %>% select(ensembl = id, symbol = Symbol, count = Blank3_count) %>% 
  mutate(sample = "GSM4584334", dataset = "GSE151536", cell_type="t_reg", level=5,
         note = "Tregs sorted from PBMC")

GSE151536 <- bind_rows(GSM4584330_Treg, GSM4584332_Treg, GSM4584334_Treg)


# Dataset GSE143690_RAW
## RCC file

# 3 Dataset GSE135390_raw_counts.csv
GSE135390_Treg17 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Treg17=lib418, BC120926_Treg17=lib429, BC121206_Treg17=lib440) %>% 
  # mutate(ensembl = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
  #                                       keys = symbol, 
  #                                       keytype = "SYMBOL", 
  #                                       column="ENSEMBL",
  #                                       multiVals = first)) %>% 
  # mutate(ensembl=unlist(ensembl)) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_reg", level=5,
         note = "CD127-CD25+CD45RA-CXCR3-CCR6+CCR4+CCR10- Treg17 sorted from peripheral blood")

GSE135390_Treg22 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Treg22=lib419, BC120926_Treg22=lib430, BC121206_Treg22=lib441) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_reg", level=5, 
         note = "CD127-CD25+CD45RA-CXCR3-CCR6+CCR4+CCR10+ Treg22 sorted from peripheral blood")

GSE135390_Treg2 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Treg2=lib420, BC120926_Treg2=lib431, BC121206_Treg2=lib442) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_reg", level=5,
         note = "CD127-CD25+CD45RA-CXCR3-CCR6-CCR4+ Treg2 sorted from peripheral blood")

GSE135390_Treg1_17 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Treg1_17=lib421, BC120926_Treg1_17=lib432, BC121206_Treg1_17=lib443) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_reg", level=5,
         note = "CD127-CD25+CD45RA-CXCR3+CCR6+ Treg1_17 sorted from peripheral blood")

GSE135390_Treg1 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Treg1=lib422, BC120926_Treg1=lib433, BC121206_Treg1=lib444) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_reg", level=5,
         note = "CD127-CD25+CD45RA-CXCR3+CCR6- Treg1 sorted from peripheral blood")


GSE135390_Th17 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th17=lib412, BC120926_Th17=lib423, BC121206_Th17=lib434) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper_h17", level=5,
         note = "CD127+CD25-CD45RA-CXCR3-CCR6+CCR4+CCR10- Th17 sorted from peripheral blood")

GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th22=lib413, BC120926_Th22=lib424, BC121206_Th22=lib435) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper_h22", level=5,
         note = "CD127+CD25-CD45RA-CXCR3-CCR6+CCR4+CCR10+ Th22 sorted from peripheral blood")

GSE135390_Th2 <- GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th2=lib414, BC120926_Th2=lib425, BC121206_Th2=lib436) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper_h2", level=5,
         note = "CD127+CD25-CD45RA-CXCR3-CCR6-CCR4+ Th2 sorted from peripheral blood")

GSE135390_raw_counts_csv %>% 
  select(symbol = X1, BC120607_Th1_17 = lib415, BC120926_Th1_17 = lib426, BC121206_Th1_17 = lib437) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper_h1_17", level=5,
         note = "CD127+CD25-CD45RA-CXCR3+CCR6+ Th1_17 sorted from peripheral blood")

GSE135390_Th1 <- GSE135390_raw_counts_csv %>% select(symbol = X1, BC120607_Th1=lib416, BC120926_Th1=lib427, BC121206_Th1=lib438) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_helper_h1", level=5,
         note = "CD127+CD25-CD45RA-CXCR3+CCR6- Th1 sorted from peripheral blood")

GSE135390_CD4_naive <- GSE135390_raw_counts_csv %>% select(symbol = X1, BC120607_Naive=lib417, BC120926_Naive=lib428, BC121206_Naive=lib439) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  filter(symbol!= "NA") %>% 
  mutate(dataset = "GSE135390", cell_type="t_CD4", level=3,
         note = "CD127+CD25-CD45RA+ naive CD4+ T cell sorted from peripheral blood")

GSE135390 <- bind_rows(GSE135390_Treg17, GSE135390_Treg22, GSE135390_Treg2, GSE135390_Treg1_17,
                       GSE135390_Treg1, GSE135390_Th17, GSE135390_Th2, GSE135390_Th1, GSE135390_CD4_naive) %>% 
  mutate(count=as.numeric(count))



# Dataset GSE148970_gene_counts_txt (NOT RAW COUNTS)
GSE148970_gene_counts_txt %>% 
  select(symbol = libID, blood_Tcm_1 = lib16824, blood_Tcm_2 = lib16834, blood_Tcm_6 = lib16882) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  filter(symbol!= "sample") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  mutate(dataset = "GSE148970", note = "CLA+CD4+CD45RA-CCR7+ T central memory sorted from peripheral blood")

GSE148970_gene_counts_txt %>% 
  select(symbol = libID, blood_Tem_1=lib16825, blood_Tem_2=lib16833, blood_Tem_3=lib16848, blood_Tem_5=lib16871, blood_Tem_6=lib16881) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  filter(symbol!= "sample") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  mutate(dataset = "GSE148970", note = "CLA+CD4+CD45RA-CCR7- T effector memory sorted from peripheral blood")

GSE148970_gene_counts_txt %>% 
  select(symbol = libID, skin_Tem_1=lib16831, skin_Tem_2=lib16840, skin_Tem_3=lib16844, skin_Tem_4=lib16856) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  filter(symbol!= "sample") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  mutate(dataset = "GSE148970", note = "CLA+CD4+CD45RA-CCR7- T effector memory sorted from processed skin")

GSE148970_gene_counts_txt %>% 
  select(symbol = libID, skin_Tcm_2=lib16841, skin_Tcm_4=lib16857) %>% 
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "count") %>% 
  filter(symbol!= "sample") %>% 
  nest(data = -sample) %>% 
  unnest(data) %>% 
  mutate(dataset = "GSE148970", note = "CLA+CD4+CD45RA-CCR7+ T central memory sorted from processed skin")


# 4 Dataset GSE138603_readcounts.txt
GSE138603_CD4_1 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_1_25_Tconv_rest_24h_GPA33_D1_ACTTGAA/genecounts.txt`) %>% 
  mutate(sample="4848_1_25", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_CD4_2 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_5_55_Tconv_rest_24h_GPA33_D6_GTGGCCT/genecounts.txt`) %>% 
  mutate(sample="4848_5_55", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_CD4_3 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_9_67_Tconv_rest_24h_metab_sort_18_ACTGATA/genecounts.txt`) %>% 
  mutate(sample="4848_9_67", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_CD4_4 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_13_79_Tconv_rest_24h_metab_sort_19_ACAGTGA/genecounts.txt`) %>% 
  mutate(sample="4848_13_79", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_CD4_5 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_2_26_Tconv_stim_24h_GPA33_D1_GATCAGA/genecounts.txt`) %>% 
  mutate(sample="4848_2_26", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")

GSE138603_CD4_6 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_6_56_Tconv_stim_24h_GPA33_D6_GTTTCGG/genecounts.txt`) %>% 
  mutate(sample="4848_6_56", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_andti-CD28_IL2")

GSE138603_CD4_7 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_10_68_Tconv_stim_24h_metab_sort_18_ATTCCTT/genecounts.txt`) %>% 
  mutate(sample="4848_10_68", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="conventional CD4+ T cell expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_andti-CD28_IL2")

GSE138603_CD4_8 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_14_80_Tconv_stim_24h_metab_sort_19_GCCAATA/genecounts.txt`) %>% 
  mutate(sample="4848_14_80", dataset="GSE138603", cell_type="t_CD4", level=3,
         note="sorted from male blood conventional CD4+ T cell expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")

GSE138603_rtreg_1 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_3_27_tTreg_rest_24h_GPA33_D1_TAGCTTA/genecounts.txt`) %>% 
  mutate(sample="4848_3_27", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_rtreg_2 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_7_57_tTreg_rest_24h_GPA33_D6_CGTACGT/genecounts.txt`) %>% 
  mutate(sample="4848_7_57", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_rtreg_3 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_11_69_tTreg_rest_24h_metab_sort_18_CGATGTA/genecounts.txt`) %>% 
  mutate(sample="4848_11_69", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_rtreg_4 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_15_81_tTreg_rest_24h_metab_sort_19_CAGATCA/genecounts.txt`) %>% 
  mutate(sample="4848_15_81", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then unstimulated for 24hrs with IL2 only")

GSE138603_atreg_5 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_4_28_tTreg_stim_24h_GPA33_D1_GGCTACA/genecounts.txt`) %>% 
  mutate(sample="4848_4_28", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")

GSE138603_atreg_6 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_8_58_tTreg_stim_24h_GPA33_D6_GAGTGGA/genecounts.txt`) %>% 
  mutate(sample="4848_8_58", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")

GSE138603_atreg_7 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_12_70_tTreg_stim_24h_metab_sort_18_TGACCAA/genecounts.txt`) %>% 
  mutate(sample="4848_12_70", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")

GSE138603_atreg_8 <- GSE138603_readcounts_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`4848_16_82_tTreg_stim_24h_metab_sort_19_CTTGTAA/genecounts.txt`) %>% 
  mutate(sample="4848_16_82", dataset="GSE138603", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs with anti-CD3_anti-CD28_IL2")


GSE138603 <- bind_rows(GSE138603_CD4_1, GSE138603_CD4_2, GSE138603_CD4_3, GSE138603_CD4_4, GSE138603_CD4_5, 
                       GSE138603_CD4_6, GSE138603_CD4_7, GSE138603_CD4_8, GSE138603_rtreg_1, GSE138603_rtreg_2,
                       GSE138603_rtreg_3, GSE138603_rtreg_4, GSE138603_atreg_5, GSE138603_atreg_6, GSE138603_atreg_7,
                       GSE138603_atreg_8)


# 5 Dataset GSE138604
GSE138604_CD4_1 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count= `5076_1_26_sort_33B_Tconv_aCD3_ATCACGA`) %>%
  mutate(sample="5076_1_26", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_CD4_2 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count= `5076_9_41_sort_29_Tconv_aCD3_GATCAGA`) %>%
  mutate(sample="5076_9_41", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_CD4_3 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count= `5076_17_56_sort_38_Tconv_aCD3_GTCCGCA`) %>%
  mutate(sample="5076_17_56", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_CD4_4 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count= `5076_2_27_sort_33B_Tconv_aCD3-28_CGATGTA`) %>%
  mutate(sample="5076_2_27", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_CD4_5 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count= `5076_10_42_sort_29_Tconv_aCD3-28_TAGCTTA`) %>%
  mutate(sample="5076_10_42", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_CD4_6 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_18_57_sort_38_Tconv_aCD3-28_GTGAAAC`) %>%
  mutate(sample="5076_18_57", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_CD4_7 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_4_29_sort_33B_Tconv_aCD3-TNFR2_TGACCAA`) %>%
  mutate(sample="5076_4_29", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")

GSE138604_CD4_8 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_12_44_sort_29_Tconv_aCD3-TNFR2_CTTGTAA`) %>%
  mutate(sample="5076_12_44", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")

GSE138604_CD4_9 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_20_59_sort_38_Tconv_aCD3-TNFR2_GTTTCGG`) %>%
  mutate(sample="5076_20_59", dataset="GSE138604", cell_type="t_CD4", level=3,
         note="sorted from male blood CD4+ conventional T cell expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")

GSE138604_aTreg_1 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_5_31_sort_33B_tTreg_aCD3_ACAGTGA`) %>%
  mutate(sample="5076_5_31", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_aTreg_2 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_13_46_sort_29_tTreg_aCD3_AGTCAAC`) %>%
  mutate(sample="5076_13_46", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_aTreg_3 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_21_61_sort_38_tTreg_aCD3_CGTACGT`) %>%
  mutate(sample="5076_21_61", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_IL2")

GSE138604_aTreg_4 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_6_32_sort_33B_tTreg_aCD3-28_GCCAATA`) %>%
  mutate(sample="5076_6_32", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_aTreg_5 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_14_47_sort_29_tTreg_aCD3-28_AGTTCCG`) %>%
  mutate(sample="5076_14_47", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_aTreg_6 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_22_62_sort_38_tTreg_aCD3-28_GAGTGGA`) %>%
  mutate(sample="5076_22_62", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-CD28_IL2")

GSE138604_aTreg_7 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_8_34_sort_33B_tTreg_aCD3-TNFR2_ACTTGAA`) %>%
  mutate(sample="5076_8_34", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")

GSE138604_aTreg_8 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_16_49_sort_29_tTreg_aCD3-TNFR2_CCGTCCC`) %>%
  mutate(sample="5076_16_49", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")

GSE138604_aTreg_9 <- GSE138604_readcounts1_txt %>% 
  select(ensembl=ensembl_gene_id, symbol=external_gene_id, count=`5076_24_64_sort_38_tTreg_aCD3-TNFR2_ATTCCTT`) %>%
  mutate(sample="5076_24_64", dataset="GSE138604", cell_type="t_reg", level=5,
         note="sorted from male blood CD4+ Tregs expanded for 1 to 2 weeks then stimulated for 24hrs by anti-CD3_anti-TNFR2_IL2")


GSE138604 <- bind_rows(GSE138604_CD4_1, GSE138604_CD4_2, GSE138604_CD4_3, GSE138604_CD4_4, GSE138604_CD4_5,
                       GSE138604_CD4_6, GSE138604_CD4_7, GSE138604_CD4_8, GSE138604_CD4_9, GSE138604_aTreg_1,
                       GSE138604_aTreg_2, GSE138604_aTreg_3, GSE138604_aTreg_4, GSE138604_aTreg_5, GSE138604_aTreg_6,
                       GSE138604_aTreg_7, GSE138604_aTreg_8, GSE138604_aTreg_9)



# Dataset GSE112770_RAW
GSM3083163 <- ReadAffy(celfile.path = "GSE112770_RAW/GSM3083163_p1427_1.CEL.gz")

# 6 Dataset GSE122941_RAW
GSM3488961_aTem <- GSM3488961_Donor_A_CSF2pos_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488961", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF+ effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSM3488962_aTem <- GSM3488962_Donor_B_CSF2pos_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488962", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF+ effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSM3488963_aTem <- GSM3488963_Donor_C_CSF2pos_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488963", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF+ effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSM3488964_aTem <- GSM3488964_Donor_A_CSF2neg_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488964", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF- effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSM3488965_aTem <- GSM3488965_Donor_B_CSF2neg_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488965", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF- effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSM3488966_aTem <- GSM3488966_Donor_C_CSF2neg_genes_counts_txt %>% 
  select(ensembl=X1, count=X2) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  mutate(sample="GSM3488966", dataset="GSE122941", cell_type="t_CD4_memory_effector", level=5,
         note="human CD4+CD25-CD45RA-CCR7-GM-CSF- effector memory T sorted from PMBC and activated by PMA_ionomycin for 3hrs")

GSE122941 <- bind_rows(GSM3488961_aTem, GSM3488962_aTem, GSM3488963_aTem, GSM3488964_aTem, GSM3488965_aTem,
                      GSM3488966_aTem)

 

# 7 Dataset GSE113891 PARSEING ERROR WHEN IMPORTING
GSE113891_Tem_1 <- GSE113891_PT_all_Count_tsv %>% 
  select(symbol=GENES, `1_2512_3_PT_TEM`, `5_2504_3_PT_TEM`,`9_2545_4_PT_TEM`, `13_1584_4_PT_TEM`, 
         `17_2488_3_PT_TEM`, `21_2490_4_PT_TEM`, `25_2496_3_PT_TEM`, `29_2511_3_PT_TEM`, `33_1858_2_PT_TEM`, 
         `37_1706_4_PT_TEM`, `41_1968_2_PT_TEM`, `45_1984_2_PT_TEM`, `49_1966_2_PT_TEM`, `53_1967_2_PT_TEM`, 
         `57_1985_2_PT_TEM`, `61_1854_2_PT_TEM`) %>% 
  pivot_longer(cols = -symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE113891", cell_type="t_CD4_memory_effector", level=5,
  note="CD3+CD4+ effector memory T cells sorted from cryopreserved PBMC of 16 donors that received vaccine boost 1-3 months 
         before experiments_Stimulation with Pertussis petides with a 24h AIM assay were used to capture and sort the antigen-specific cells")

GSE113891_Tem_2 <- GSE113891_PT_all_Count_tsv %>% 
  select(symbol=GENES, `3_2512_3_TT_TEM`, `7_2504_3_TT_TEM`, `11_2545_4_TT_TEM`, `15_1584_4_TT_TEM`, `19_2488_3_TT_TEM`, 
         `23_2490_4_TT_TEM`, `27_2496_3_TT_TEM`, `31_2511_3_TT_TEM`, `35_1858_2_TT_TEM`, `39_1706_4_TT_TEM`, 
         `43_1968_2_TT_TEM`, `47_1984_2_TT_TEM`, `51_1966_2_TT_TEM`, `55_1967_2_TT_TEM`, `59_1985_2_TT_TEM`, 
         `63_1854_2_TT_TEM`) %>% 
  pivot_longer(cols = -symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE113891", cell_type="t_CD4_memory_effector", level=5,
  note="CD3+CD4+ effector memory T cells sorted from cryopreserved PBMC of 16 donors that received vaccine boost 1-3 months 
         before experiments_Stimulation with Tetanus petides with a 24h AIM assay were used to capture and sort the antigen-specific cells")

GSE113891_Tcm_1 <- GSE113891_PT_all_Count_tsv %>% 
  select(symbol=GENES, `2_2512_3_PT_TCM`, `6_2504_3_PT_TCM`, `10_2545_4_PT_TCM`, `14_1584_4_PT_TCM`, `18_2488_3_PT_TCM`, 
         `22_2490_4_PT_TCM`, `26_2496_3_PT_TCM`, `30_2511_3_PT_TCM`, `34_1858_2_PT_TCM`, `38_1706_4_PT_TCM`, 
         `42_1968_2_PT_TCM`, `46_1984_2_PT_TCM`, `50_1966_2_PT_TCM`, `54_1967_2_PT_TCM`, `58_1985_2_PT_TCM`, 
         `62_1854_2_PT_TCM`) %>% 
  pivot_longer(cols = -symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE113891", cell_type="t_CD4_memory_central", level=5, 
  note="CD3+CD4+ central memory T cells sorted from cryopreserved PBMC of 16 donors that received vaccine boost 1-3 months 
         before experiments_Stimulation with Pertussis petides with a 24h AIM assay were used to capture and sort the antigen-specific cells")

GSE113891_Tcm_2 <- GSE113891_PT_all_Count_tsv %>% 
  select(symbol=GENES, `4_2512_3_TT_TCM`, `8_2504_3_TT_TCM`, `12_2545_4_TT_TCM`, `16_1584_4_TT_TCM`, `20_2488_3_TT_TCM`, 
         `24_2490_4_TT_TCM`, `28_2496_3_TT_TCM`, `32_2511_3_TT_TCM`, `36_1858_2_TT_TCM`, `40_1706_4_TT_TCM`, 
         `44_1968_2_TT_TCM`, `48_1984_2_TT_TCM`, `52_1966_2_TT_TCM`, `56_1967_2_TT_TCM`, `60_1985_2_TT_TCM`, 
         `64_1854_2_TT_TCM`) %>% 
  pivot_longer(cols = -symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE113891", cell_type="t_CD4_memory_central", level=5,
  note="CD3+CD4+ central memory T cells sorted from cryopreserved PBMC of 16 donors that received vaccine boost 1-3 months 
         before experiments_Stimulation with Tetanus petides with a 24h AIM assay were used to capture and sort the antigen-specific cells")

GSE113891 <- bind_rows(GSE113891_Tem_1, GSE113891_Tem_2, GSE113891_Tcm_1, GSE113891_Tcm_2)



# 8 Dataset GSE85294_read_count_table.txt
GSE85294_naive_CD4 <- GSE85294_read_count_table_txt %>% 
  separate(Feature, c("ensembl", "symbol"), sep="_", extra="merge") %>% 
  select(ensembl, symbol, GSM2264211=`3026-TA-1`, GSM2264214=`3026-TA-4`, GSM2264217=`3026-TA-7`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE85294", cell_type="t_CD4_naive", level=3,
         note="CD45RA_hi CCR7_hi naive CD4+ T cells sorted from PMBC")

GSE85294_Tcm <- GSE85294_read_count_table_txt %>% 
  separate(Feature, c("ensembl", "symbol"), sep="_", extra="merge") %>% 
  select(ensembl, symbol, GSM2264212=`3026-TA-2`, GSM2264215=`3026-TA-5`, GSM2264218=`3026-TA-8`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE85294", cell_type="t_CD4_memory_central", level=5,
         note="CD45RA_lo CCR7_hi CD4+ central memory T cells sorted from PMBC")

GSE85294_Tem <- GSE85294_read_count_table_txt %>% 
  separate(Feature, c("ensembl", "symbol"), sep="_", extra="merge") %>% 
  select(ensembl, symbol, GSM2264213=`3026-TA-3`, GSM2264216=`3026-TA-6`, GSM2264219=`3026-TA-9`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE85294", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RA_lo CCR7_lo CD4+ effector memory T cells sorted from PMBC")


GSE85294 <- bind_rows(GSE85294_naive_CD4, GSE85294_Tcm, GSE85294_Tem)



# 9 Dataset GSE89404_RAW
GSM2370618_naive_CD4 <- GSM2370618_CD4_T_cells_Naive_0h_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370618", dataset="GSE89404", cell_type="t_CD4_naive", level=3,
         note="CD45RO-CD27+ naive resting CD4+ T cells sorted from peripheral blood ")

GSM2370619_naive_CD4 <- GSM2370619_CD4_T_cells_Naive_40m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370619", dataset="GSE89404", cell_type="t_CD4_naive", level=3,
         note="CD45RO-CD27+ naive CD4+ T cells from PB activated for 40mins by anti-CD3/28 beads")

GSM2370620_naive_CD4 <- GSM2370620_CD4_T_cells_Naive_150m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370620", dataset="GSE89404", cell_type="t_CD4_naive", level=3,
         note="CD45RO-CD27+ naive CD4+ T cells from PB activated for 150mins by anti-CD3/28 beads")

GSM2370621_naive_CD4 <- GSM2370621_CD4_T_cells_Naive_15h_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370621", dataset="GSE89404", cell_type="t_CD4_memory_central", level=5,
         note="CD45RO-CD27+ naive CD4+ T cells from PB activated for 15hrs by anti-CD3/28 beads")

GSM2370622_Tcm <- GSM2370622_CD4_T_cells_CM_0h_csv %>%  select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370622", dataset="GSE89404", cell_type="t_CD4_memory_central", level=5,
         note="CD45RO+CD27+ resting CD4+ central memory T cells sorted from peripheral blood")

GSM2370623_Tcm <- GSM2370623_CD4_T_cells_CM_40m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370623", dataset="GSE89404", cell_type="t_CD4_memory_central", level=5,
         note="CD45RO+CD27+ CD4+ central memory T cells from PB activated for 40mins by anti-CD3/28 beads")

GSM2370624_Tcm <- GSM2370624_CD4_T_cells_CM_150m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370624", dataset="GSE89404", cell_type="t_CD4_memory_central", level=5,
         note="CD45RO+CD27+ CD4+ central memory T cells from PB activated for 150mins by anti-CD3/28 beads")

GSM2370625_Tcm <- GSM2370625_CD4_T_cells_CM_15h_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370625", dataset="GSE89404", cell_type="t_CD4_memory_central", level=5,
         note="CD45RO+CD27+ CD4+ central memory T cells from PB activated for 15hrs by anti-CD3/28 beads")

GSM2370626_Tem <- GSM2370626_CD4_T_cells_EM_0h_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370626", dataset="GSE89404", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RO+CD27- resting CD4+ effector memory T cells sorted from peripheral blood")

GSM2370627_Tem <- GSM2370627_CD4_T_cells_EM_40m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370627", dataset="GSE89404", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RO+CD27- CD4+ effector memory T cells from PB activated for 40mins by anti-CD3/28 beads")

GSM2370628_Tem <- GSM2370628_CD4_T_cells_EM_150m_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370628", dataset="GSE89404", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RO+CD27- CD4+ effector memory T cells from PB activated for 150mins by anti-CD3/28 beads")

GSM2370629_Tem <- GSM2370629_CD4_T_cells_EM_15h_csv %>% select(ensembl=refseq_id, symbol=gene_id, count=TOT_R_0) %>% 
  mutate(sample="GSM2370629", dataset="GSE89404", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RO+CD27- CD4+ effector memory T cells from PB activated for 15hrs by anti-CD3/28 beads")

GSE89404 <- bind_rows(GSM2370618_naive_CD4, GSM2370619_naive_CD4, GSM2370620_naive_CD4, GSM2370621_naive_CD4,
                      GSM2370622_Tcm, GSM2370623_Tcm, GSM2370624_Tcm, GSM2370625_Tcm, GSM2370626_Tem ,
                      GSM2370627_Tem, GSM2370628_Tem, GSM2370629_Tem)



# Mast cells

# 10 Dataset GSE125887_RAW  MORE THAN ONE COLUMN OF COUNTS
GSM3583947_mast <- GSM3583947_RHM5399_star_hg19_untreated_rep1_counts_tsv %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  slice(-c(1, 2, 3)) %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583947", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="untreated peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

View(GSM3583947_RHM5399_star_hg19_untreated_rep1_counts_tsv)

GSM3583950_mast <- GSM3583950_RHM5402_star_hg19_untreated_rep2_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583950", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="untreated peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583953_mast <- GSM3583953_RHM5405_star_hg19_untreated_rep3_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583953", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="untreated peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583956_mast <- GSM3583956_RHM5408_star_hg19_untreated_rep4_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583953", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="untreated peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583948_mast <- GSM3583948_RHM5400_star_hg19_IgE_rep1_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583948", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583951_mast <- GSM3583951_RHM5403_star_hg19_IgE_rep2_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583951", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583954_mast <- GSM3583954_RHM5406_star_hg19_IgE_rep3_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583954", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583957_mast <- GSM3583957_RHM5409_star_hg19_IgE_rep4_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583957", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583949_mast <- GSM3583949_RHM5401_star_hg19_anti_IgE_rep1_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583949", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="2hr anti-IgE activation of overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583952_mast <- GSM3583952_RHM5404_star_hg19_anti_IgE_rep2_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583952", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="2hr anti-IgE activation of overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583955_mast <- GSM3583955_RHM5407_star_hg19_anti_IgE_rep3_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583955", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="2hr anti-IgE activation of overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")

GSM3583958_mast <- GSM3583958_RHM5410_star_hg19_anti_IgE_rep4_counts_tsv %>% 
  slice(-c(1, 2, 3)) %>% 
  mutate(symbol=N_unmapped, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sequence_run", values_to="count") %>% 
  mutate(sample="GSM3583958", dataset="GSE125887", cell_type="mast_cell", level=2,
         note="2hr anti-IgE activation of overnight IgE sensitised peripheral blood-derived mast cells (CD34+ cells from buffy coats cultured for 10 weeks in hSCF_hIL-6 hIL-3(first 3 weeks only))")


GSE125887 <- bind_rows(GSM3583947_mast, GSM3583950_mast, GSM3583953_mast, GSM3583956_mast, GSM3583948_mast, 
                       GSM3583951_mast, GSM3583954_mast, GSM3583957_mast, GSM3583949_mast, GSM3583952_mast,
                       GSM3583955_mast, GSM3583958_mast)





# Dataset GSE71247_RAW # NOT RAW COUNT
GSM1831353_mast <- GSM1831353_072415_PB_MC_PC_txt %>% select(symbol=`Gene Symbol`, count=Value) %>% 
  mutate(sample="GSM3583956", dataset="GSE71247", cell_type="mast_cell", level=2,
         note="PB-derived cultured mast cells generated from CD34+ hematopoietic stem cells in the PB with IL-3 IL-6 and stem cell factor")

View(GSM1831353_mast)

# 11 Dataset GSE115103_raw_counts
GSE115103_Th1 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165169=S3855, GSM3165179=S3868, GSM3165189=S3882) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper_h1", level=5,
         note="CD4+ Th1 cells sorted from PBMC")

GSE115103_Th2 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165170=S3856, GSM3165180=S3869, GSM3165190=S3883) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper_h2", level=5,
         note="CD4+ Th2 cells sorted from PBMC")

GSE115103_Th17 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165171=S3857, GSM3165181=S3870, GSM3165191=S3884) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper_h17", level=5,
         note="CD4+ Th17 cells sorted from PBMC")


GSE115103_Th1_17 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165172=S3858, GSM3165182=S3871, GSM3165192=S3885) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper_h1_17", level=5,
         note="CD4+ Th1_Th17 cells sorted from PBMC")

GSE115103_Th22 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165173=S3859, GSM3165183=S3872, GSM3165193=S3886) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_helper_h22", level=5,
         note="CD4+ Th22 cells sorted from PBMC")

GSE115103_Tc1 <- GSE115103_raw_counts_txt %>%
  select(ensembl=X1, GSM3165164=S3848, GSM3165174=S3861, GSM3165184=S3875) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_CD8", level=3,
         note="CD8+ Tc1 cells sorted from PBMC")

GSE115103_Tc2 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165165=S3849, GSM3165175=S3862, GSM3165185=S3876) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_CD8", level=3,
         note="CD8+ Tc2 cells sorted from PBMC")

GSE115103_Tc17 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165166=S3850, GSM3165176=S3863, GSM3165186=S3877) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_CD8", level=3,
         note="CD8+ Tc17 cells sorted from PBMC")

GSE115103_Tc1_17 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165167=S3851, GSM3165177=S3864, GSM3165187=S3878) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_CD8", level=3,
         note="CD8+ Tc1_Tc17 cells sorted from PBMC")

GSE115103_Tc22 <- GSE115103_raw_counts_txt %>% 
  select(ensembl=X1, GSM3165168=S3852, GSM3165178=S3865, GSM3165188=S3879) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>%  
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115103", cell_type="t_CD8", level=3,
         note="CD8+ Tc22 cells sorted from PBMC")


GSE115103 <- bind_rows(GSE115103_Th1, GSE115103_Th2, GSE115103_Tc1, GSE115103_Tc2, GSE115103_Tc17, GSE115103_Tc1_17,
                       GSE115103_Tc22)





# 12 Dataset GSE151079_RNAseq_thymocytes_raw_counts
GSE151079_Tgd_CD1pos <- GSE151079_RNAseq_thymocytes_raw_counts_txt %>% 
  select(ensembl=X1, GSM4566068=`gd CD1+ R1`, GSM4566079=`gd CD1+ R2`) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE151079", cell_type="t_gamma_delta", level=3,
         note="TCRgd+ CD3+ CD1+ thymocytes from post natal thymus")

GSE151079_Tgd_CD1neg <- GSE151079_RNAseq_thymocytes_raw_counts_txt %>% 
  select(ensembl=X1, symbol=symbol, GSM4566068=`gd CD1- R1`, GSM4566079=`gd CD1- R2` ) %>% 
  ensembl_to_symbol(ensembl) %>%
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols = -c(ensembl, symbol), names_to="sample", values_to="count") %>%
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE151079", cell_type="t_gamma_delta", level=3,
         note="TCRgd+ CD3+ CD1- thymocytes from post natal thymus")


GSE151079 <- bind_rows(GSE151079_Tgd_CD1pos, GSE151079_Tgd_CD1neg)





# 13  Dataset GSE107981_hu_counts
GSE107981_CD4_naive <- GSE107981_hu_counts_txt %>% 
  select(symbol=Gene_ID, GSM2885044=Mahu001, GSM2885045=Mahu002, GSM2885046=Mahu003, GSM2885047=Mahu004) %>% 
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE107981", cell_type="t_CD4", level=3,
         note="CD45RA hi CD45RO- naive CD4+ freshly isolated from PBMC")

GSE107981_Th1 <- GSE107981_hu_counts_txt %>% 
  select(symbol=Gene_ID, GSM2885048=Mahu005, GSM2885049=Mahu006, GSM2885050=Mahu007, GSM2885051=Mahu008) %>% 
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE107981", cell_type="t_CD4", level=3,
         note="CD45RA hi CD45RO- CXCR3+CCR4-CCR6- CD4+ Th1 freshly isolated from PBMC")


GSE107981 <- bind_rows(GSE107981_CD4_naive, GSE107981_Th1)



# 14 Dataset GSE75011_Raw_counts
GSE75011_Th2 <- GSE75011_Raw_counts_tsv %>% 
  select(symbol=X1, GSM1940714=D66_HC, GSM1940715=D67_HC, GSM1940716=D68_HC, GSM1940717=D69_HC, GSM1940718=D70_HC, 
         GSM1940719=D71_HC, GSM1940720=D72_HC, GSM1940721=D73_HC, GSM1940722=D74_HC, GSM1940723=D75_HC, 
         GSM1940724=D76_HC, GSM1940725=D77_HC, GSM1940726=D78_HC, GSM1940727=D79_HC, GSM1940728=D80_HC) %>% 
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE75011", cell_type="t_helper_h2", level=5,
         note="CD3+CD4+CD45RA-CCR7+CD25-CCR4+ resting Th2 cells freshly sorted from PBMCs")

GSE75011 <- GSE75011_Th2




# 15  Dataset GSE118974_RAW
GSM3355220_Th17 <- GSM3355220_180036_MK198_6_Th17_D1_S85_L005_txt %>% 
  select(symbol=`#`, count=`"-g"`) %>% slice(-1) %>% 
  mutate(sample="GSM3355220", dataset="GSE118974", cell_type="t_helper_h17", level=5,
         note="naive CD4+ PBMCs from umbilical cord blood of healthy neonates cultured for 72hrs under Th17 polarising condition IL-6 + IL1-β + TGF-β + anti-IFNγ + anti-IL4 + anti-CD3 + anti-CD28")

GSM3355221_Th17 <- GSM3355221_180036_MK198_6_Th17_D2_S84_L005_txt %>% 
  select(symbol=`#`, count=`"-g"`) %>% slice(-1) %>% 
  mutate(sample="GSM3355221", dataset="GSE118974", cell_type="t_helper_h17", level=5,
         note="naive CD4+ PBMCs from umbilical cord blood of healthy neonates cultured for 72hrs under Th17 polarising condition IL-6 + IL1-β + TGF-β + anti-IFNγ + anti-IL4 + anti-CD3 + anti-CD28")

GSM3355222_Th17 <- GSM3355222_180036_MK198_6_Th17_D3_S83_L005_txt %>% 
  select(symbol=`#`, count=`"-g"`) %>% slice(-1) %>% 
  mutate(sample="GSM3355222", dataset="GSE118974", cell_type="t_helper_h17", level=5,
         note="naive CD4+ PBMCs from umbilical cord blood of healthy neonates cultured for 72hrs under Th17 polarising condition IL-6 + IL1-β + TGF-β + anti-IFNγ + anti-IL4 + anti-CD3 + anti-CD28")

GSM3355223_Th17 <- GSM3355223_180036_MK198_6_Th17_D4_S82_L005_txt %>% 
  select(symbol=`#`, count=`"-g"`) %>% slice(-1) %>% 
  mutate(sample="GSM3355223", dataset="GSE118974", cell_type="t_helper_h17", level=5,
         note="naive CD4+ PBMCs from umbilical cord blood of healthy neonates cultured for 72hrs under Th17 polarising condition IL-6 + IL1-β + TGF-β + anti-IFNγ + anti-IL4 + anti-CD3 + anti-CD28")

GSM3355224_Th17 <- GSM3355224_180036_MK198_6_Th17_D5_S87_L005_txt %>%
  select(symbol=`#`, count=`"-g"`) %>% slice(-1) %>% 
  mutate(sample="GSM3355224", dataset="GSE118974", cell_type="t_helper_h17", level=5,
         note="naive CD4+ PBMCs from umbilical cord blood of healthy neonates cultured for 72hrs under Th17 polarising condition IL-6 + IL1-β + TGF-β + anti-IFNγ + anti-IL4 + anti-CD3 + anti-CD28")


GSE118974 <- bind_rows(GSM3355220_Th17, GSM3355221_Th17, GSM3355222_Th17, GSM3355223_Th17, GSM3355224_Th17) %>% 
  mutate(count=as.numeric(count))

 


# 16 Dataset GSE155715 ENTREZ ID INSTEAD OF ENSEMBLE ID!
GSE155715_aCD8 <- GSE155715_Raw_gene_counts_matrix_txt %>% 
  select(ensembl=Geneid, GSM4711152=`Donor A_vehicle`, GSM4711154=`Donor B_vehicle`, GSM4711156=`Donor C_vehicle`) %>% 
  mutate(ensembl=as.character(ensembl)) %>% 
  mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                        keys = ensembl, 
                                        keytype = "ENTREZID", 
                                        column="SYMBOL", 
                                        multiVals = "first")) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE155715", cell_type="t_CD8", level=3,
         note="naive CD8+ T cells activated by aCD3 beads and anti-CD28 for 4 days then cultured in IL-2 supplied media for another 8 days")

GSE155715 <- GSE155715_aCD8



# 17 Dataset GSE147394
GSE147394_Tcm <- GSE147394_datamatrix_rawcounts_dataset1_txt %>% 
  select(ensembl=gene_id, GSM4429943=HD1_TCM, GSM4429946=HD2_TCM, GSM4429949=HD3_TCM, 
         GSM4429953=HD4_TCM, GSM4429956=HD5_TCM) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[1-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE147394", cell_type="t_CD8_memory_central", level=5,
         note="Tigit-PD1- CD8+ T central memory cells sorted from PBMCs")

GSE147394_Tem <- GSE147394_datamatrix_rawcounts_dataset1_txt %>% 
  select(ensembl=gene_id, GSM4429947=HD2_TEM, GSM4429951=HD3_TEM, GSM4429958=HD5_TEM) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[1-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE147394", cell_type="t_CD8_memory_effector", level=5,
         note="CD8+ T effector memory cells sorted from PBMCs")


GSE147394 <- bind_rows(GSE147394_Tcm, GSE147394_Tem)





# 18 Dataset GSE136200
GSE136200_CD8 <- GSE136200_CD8_linear_raw_txt %>% 
  select(symbol=Gene, GSM4042527=H43, GSM4042528=H44, GSM4042529=H45, GSM4042530=H46) %>% 
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE136200", cell_type="t_CD8", level=3,
         note="CD8 from PBMC")

GSE136200 <- GSE136200_CD8 



# 19 Dataset GSE118829_RAW
GSE118829_CD4TCM_txt <- GSM3348321_HC1_CD4TCM_txt %>% 
  mutate(symbol=Gene, GSM3348321=Count, .keep="unused") %>% 
  mutate(GSM3348328=GSM3348328_HC10_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348335=GSM3348335_HC2_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348341=GSM3348341_HC3_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348348=GSM3348348_HC4_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348355=GSM3348355_HC5_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348362=GSM3348362_HC6_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348369=GSM3348369_HC7_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348376=GSM3348376_HC8_CD4TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348383=GSM3348383_HC9_CD4TCM_txt$Count, .keep="unused") %>% 
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD4_memory_central", level=5,
         note="CD45RA-CCR7+ CD4+ central memory T cells from PBMC")

GSE118829_CD4TEM_txt <- GSM3348322_HC1_CD4TEM_txt %>% 
  mutate(symbol=Gene, GSM3348322=Count, .keep="unused") %>% 
  mutate(GSM3348329=GSM3348329_HC10_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348336=GSM3348336_HC2_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348342=GSM3348342_HC3_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348349=GSM3348349_HC4_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348356=GSM3348356_HC5_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348363=GSM3348363_HC6_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348370=GSM3348370_HC7_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348377=GSM3348377_HC8_CD4TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348384=GSM3348384_HC9_CD4TEM_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD4_memory_effector", level=5,
         note="CD45RA-CCR7- CD4+ effector memory T cells from PBMC")

GSE118829_CD4_naive_txt <- GSM3348323_HC1_CD4TN_txt %>% 
  mutate(symbol=Gene, GSM3348323=Count, .keep="unused") %>% 
  mutate(GSM3348330=GSM3348330_HC10_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348337=GSM3348337_HC2_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348343=GSM3348343_HC3_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348350=GSM3348350_HC4_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348357=GSM3348357_HC5_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348364=GSM3348364_HC6_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348371=GSM3348371_HC7_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348378=GSM3348378_HC8_CD4TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348385=GSM3348385_HC9_CD4TN_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD4", level=3,
         note="CD45RA+CCR7+ CD4+ naive T cells from PBMC")

GSE118829_CD8_TCM_txt <- GSM3348324_HC1_CD8TCM_txt %>% 
  mutate(symbol=Gene, GSM3348324=Count, .keep="unused") %>% 
  mutate(GSM3348324=GSM3348324_HC1_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348331=GSM3348331_HC10_CD8TCM_txt$Count, .keep="unused") %>% 
  mutate(GSM3348338=GSM3348338_HC2_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348344=GSM3348344_HC3_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348351=GSM3348351_HC4_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348358=GSM3348358_HC5_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348365=GSM3348365_HC6_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348372=GSM3348372_HC7_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348379=GSM3348379_HC8_CD8TCM_txt$Count, .keep="unused") %>%
  mutate(GSM3348386=GSM3348386_HC9_CD8TCM_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD8_memory_central", level=5,
         note="CD45RA-CCR7+ CD8+ central memory cells from PBMC")
  
GSE118829_CD8_TEM_txt <- GSM3348325_HC1_CD8TEM_txt %>% 
  mutate(symbol=Gene, GSM3348325=Count, .keep="unused") %>% 
  mutate(GSM3348332=GSM3348332_HC10_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348339=GSM3348339_HC2_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348345=GSM3348345_HC3_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348352=GSM3348352_HC4_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348359=GSM3348359_HC5_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348366=GSM3348366_HC6_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348373=GSM3348373_HC7_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348380=GSM3348380_HC8_CD8TEM_txt$Count, .keep="unused") %>%
  mutate(GSM3348387=GSM3348387_HC9_CD8TEM_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD8_memory_effector", level=5,
         note="CD45RA-CCR7- CD8+ effector memory cells from PBMC")
  
GSE118829_CD8_naive_txt <- GSM3348327_HC1_CD8TN_txt %>% 
  mutate(symbol=Gene, GSM3348327=Count, .keep="unused") %>% 
  mutate(GSM3348334=GSM3348334_HC10_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348340=GSM3348340_HC2_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348347=GSM3348347_HC3_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348354=GSM3348354_HC4_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348361=GSM3348361_HC5_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348368=GSM3348368_HC6_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348375=GSM3348375_HC7_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348382=GSM3348382_HC8_CD8TN_txt$Count, .keep="unused") %>%
  mutate(GSM3348389=GSM3348389_HC9_CD8TN_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD8", level=3,
         note="CD45RA+CCR7+ naive CD8+ T cells from PBMC")
  
  
GSE118829_CD8_TEMRA_txt <- GSM3348326_HC1_CD8TEMRA_txt %>% 
  mutate(symbol=Gene, GSM3348326=Count, .keep="unused") %>% 
  mutate(GSM3348333=GSM3348333_HC10_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348346=GSM3348346_HC3_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348353=GSM3348353_HC4_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348360=GSM3348360_HC5_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348367=GSM3348367_HC6_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348374=GSM3348374_HC7_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348381=GSM3348381_HC8_CD8TEMRA_txt$Count, .keep="unused") %>%
  mutate(GSM3348388=GSM3348388_HC9_CD8TEMRA_txt$Count, .keep="unused") %>%
  pivot_longer(cols=-symbol, names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE118829", cell_type="t_CD8_memory_effector", level=5,
         note="CD45RA+CCR7- CD8+ CD45RA+ T effector memory cells from PBMC")


GSE118829 <- bind_rows(GSE118829_CD4TCM_txt, GSE118829_CD4TEM_txt, GSE118829_CD4_naive_txt, 
                       GSE118829_CD8_TCM_txt, GSE118829_CD8_TEM_txt, GSE118829_CD8_naive_txt, 
                       GSE118829_CD8_TEMRA_txt)




# Dataset GSE134576_RAW
GSE134576_CD8_TEMRA_txt <- untar("GSE134576_RAW/GSM3956368_T3.tar.gz")
rm(GSE134576_CD8_TEMRA_txt)  



# 20 Dataset GSE123812_CD4_bulk_RNA_counts
GSE123812_CD4_naive <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, GSM3511697=Naive_1, GSM3511698=Naive_2, GSM3511699=Naive_3, GSM3511700=Naive_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_CD4", level=3,
         note="CD25-CD45RA+ naive CD4+ primary T cells sorted from PBMCs")

GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, GSM3511701=`Tfh1-17_1`, GSM3511702=`Tfh1-17_2`, GSM3511703=`Tfh1-17_3`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_follicular_helper_h1_17", level=5,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3+CCR6+ CD4+ primary follicular Th1-17 cells sorted from PBMCs")

GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511704=Tfh1_1, GSM3511705=Tfh1_2, GSM3511706=Tfh1_3, GSM3511707=Tfh1_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_follicular_helper_h1", level=5,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3+CCR6- CD4+ follicular Th1 cells sorted from PBMCs")

GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511708=Tfh17_1, GSM3511709=Tfh17_2, GSM3511710=Tfh17_3, GSM3511711=Tfh17_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_follicular_helper_17", level=5,
         note="CXCR5+ CD25-IL7RhiCD45RA-CXCR3-CCR6+ CD4+ follicular Th17 cells sorted from PBMCs")

GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511712=Tfh2_1, GSM3511713=Tfh2_2, GSM3511714=Tfh2_3, GSM3511715=Tfh2_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_follicular_helper_h2", level=5,
         note="CXCR5+ CD4+CD25-IL7RhiCD45RA-CXCR3-CCR6- CD4+ follicular Th2 cells sorted from PBMCs")

GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511716=`Th1-17_1`, GSM3511717=`Th1-17_2`, GSM3511718=`Th1-17_3`, GSM3511719=`Th1-17_4`) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper_h1_17", level=5,
         note="CD25-IL7RhiCD45RA-CXCR3+CCR6+ CD4+ Th1-17 cells sorted from PBMCs")

GSE123812_Th1 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511720=Th1_1, GSM3511721=Th1_2, GSM3511722=Th1_3, GSM3511723=Th1_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper_h1", level=5,
         note="CD25-IL7RhiCD45RA-CXCR3+CCR6- CD4+ Th1 cells sorted from PBMCs")

GSE123812_Th17 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511724=Th17_1, GSM3511725=Th17_2, GSM3511726=Th17_3, GSM3511727=Th17_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper_h17", level=5,
         note="CD25-IL7RhiCD45RA-CXCR3-CCR6+ CD4+ Th17 cells sorted from PBMCs")

GSE123812_Th2 <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511728=Th2_1, GSM3511729=Th2_2, GSM3511730=Th2_3, GSM3511731=Th2_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_helper_h2", level=5,
         note=" CD4+ CD25-IL7RhiCD45RA-CXCR3-CCR6- Th2 cells sorted from PBMCs")

GSE123812_rTreg <- GSE123812_CD4_bulk_RNA_counts_txt %>%  
  select(ensembl=EnsemblGeneID, symbol=GeneName, 
         GSM3511732=Treg_1, GSM3511733=Treg_3, GSM3511734=Treg_4) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE123812", cell_type="t_reg", level=5,
         note=" CD4+ CD25+IL7Rlo Treg cells sorted from PBMCs")



GSE123812 <- bind_rows(GSE123812_CD4_naive, GSE123812_Th1, GSE123812_Th17, GSE123812_Th2, GSE123812_rTreg)




# 21 Dataset GSE153104_all_counts
GSE153104_CD4_memory <- GSE153104_all_counts_tsv %>% 
  select(ensembl=gene_id, GSM4634224=`2950_CD4`, GSM4634225=`2951_CD4`, GSM4634226=`2965_CD4`, GSM4634227=`2964_CD4`,
         GSM4634228=`2999_CD4`, GSM4634229=`3015_CD4`, GSM4634230=`2972_CD4`, GSM4634238=`3169_CD4`, GSM4634244=`3292_CD4`,
         GSM4634245=`3344_CD4`, GSM4634248=`3417_CD4`, GSM4634249=`3454_CD4`, GSM4634250=`3420_CD4`, GSM4634251=`3635_CD4`,
         GSM4634252=`3677_CD4`, GSM4634253=`2930_CD4`, GSM4634254=`2948_CD4`, GSM4634255=`2978_CD4`, GSM4634260=`3831_CD4`, 
         GSM4634270=`3180_CD4m`, GSM4634271=`3181_CD4m`, GSM4634272=`3190_CD4m`, GSM4634273=`3195_CD4m`, GSM4634274=`3204_CD4m`,
         GSM4634275=`3207_CD4m`, GSM4634276=`3277_CD4m`, GSM4634277=`3278_CD4m`, GSM4634278=`3345_CD4m`) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[1-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>%
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE153104", cell_type="t_CD4_memory", level=4,
         note="CD4+ memory T cells are stained and sorted from thawed frozen PBMCs")

GSE153104_CD8_memory <- GSE153104_all_counts_tsv %>% 
  select(ensembl=gene_id, GSM4634279=`2950_CD8`, GSM4634280=`2951_CD8`, GSM4634281=`2965_CD8`, GSM4634282=`2964_CD8`, 
         GSM4634283=`2999_CD8`, GSM4634284=`3015_CD8`, GSM4634285=`2972_CD8`, GSM4634286=`3019_CD8`, GSM4634293=`3169_CD8`,
         GSM4634299=`3292_CD8`, GSM4634300=`3344_CD8`, GSM4634303=`3417_CD8`, GSM4634304=`3454_CD8`, GSM4634305=`3420_CD8`,
         GSM4634306=`3635_CD8`, GSM4634307=`3677_CD8`, GSM4634308=`2930_CD8`, GSM4634309=`2948_CD8`, GSM4634310=`2978_CD8`,
         GSM4634315=`3831_CD8`, GSM4634325=`3180_CD8m`, GSM4634326=`3181_CD8m`, GSM4634327=`3190_CD8m`, GSM4634328=`3195_CD8m`,
         GSM4634329=`3197_CD8m`, GSM4634330=`3204_CD8m`, GSM4634331=`3207_CD8m`, GSM4634332=`3277_CD8m`, GSM4634333=`3278_CD8m`,
         GSM4634334=`3345_CD8m`) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[1-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>%
  mutate(symbol=transcript, .keep="unused") %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE153104", cell_type="t_CD8_memory", level=4,
         note="CD8+ memory T cells are stained and sorted from thawed frozen PBMCs")

GSE153104 <- bind_rows(GSE153104_CD4_memory, GSE153104_CD8_memory)




# 22 Dataset GSE133527 ENTREZ ID!!!
GSE133527_M2_Macrophage <- GSE133527_gene_counts_1_tsv %>% 
  select(ensembl=X1, GSM3911300=ctrl_81, GSM3911301=ctrl_85, GSM3911302=ctrl_95, GSM3911303=ctrl_50, 
         GSM3911304=ctrl_52) %>% 
  mutate(ensembl=as.character(ensembl)) %>%
  mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                        keys = as.character(ensembl), 
                                        keytype = "ENTREZID", 
                                        column="SYMBOL", 
                                        multiVals = "first")) %>% 
  pivot_longer(cols=-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE133527", cell_type="macrophage_M2", level=4,
         note="in vitro derived M2 macrophages")
  
GSE133527 <- GSE133527_M2_Macrophage




# 23 GSE131792_RawCounts.csv
GSE131792_141hi_cDC <- GSE131792_RawCounts_csv %>% 
  select(ensembl=X1, GSM3819865=`NBM-11-CD141Hi`, GSM3819868=`NBM-2-CD141Hi`, GSM3819871=`NBM-4-CD141Hi`, 
         GSM3819874=`NBM-6-CD141Hi`, GSM3819877=`NBM-7-CD141Hi`, GSM3819880=`NBM-9-CD141Hi`) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE131792", cell_type="dendritic_myeloid_immature", level=4,
         note="CD141hi bone marrow cells")

GSE131792_CD1c_cDC <- GSE131792_RawCounts_csv %>% 
  select(ensembl=X1, GSM3819866=`NBM-11-CD1c`, GSM3819869=`NBM-2-CD1c`, GSM3819872=`NBM-4-CD1c`, 
         GSM3819875=`NBM-6-CD1c`, GSM3819878=`NBM-7-CD1c`, GSM3819881=`NBM-9-CD1c`) %>% 
  extract(ensembl, "ensembl", "([a-zA-Z0-9]+)\\.[0-9]*") %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE131792", cell_type="dendritic_myeloid_immature", level=4,
         note="CD1c+ bone marrow cells")


GSE131792 <- bind_rows(GSE131792_141hi_cDC, GSE131792_CD1c_cDC)



# 24 Dataset GSE115736_Haemopedia-Human-RNASeq_raw.txt
GSE115736_CD123myDC <- GSE115736_Haemopedia_Human_RNASeq_raw_txt %>% 
  select(ensembl=X1, GSM3188494=myDC123.1, GSM3188495=myDC123.2) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115736", cell_type="dendritic_myeloid_immature", level=4,
         note="CD3-CD14-BDCA2-BDCA1+CD123+ myeloid DC from PBMC")

GSE115736_myDC <- GSE115736_Haemopedia_Human_RNASeq_raw_txt %>% 
  select(ensembl=X1, GSM3188505=myDC.1, GSM3188515=myDC.2, GSM3188534=myDC.3) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE115736", cell_type="dendritic_myeloid_immature", level=4,
         note="CD3-CD14-BDCA2-BDCA1+CD123- myeloid DC from PBMC")


GSE115736 <- bind_rows(GSE115736_CD123myDC, GSE115736_myDC)





# 25 Dataset GSE70106_counts.txt
GSE70106_mDC <- GSE70106_counts_txt %>% 
  select(ensembl=Feature, GSM1717152=`SM_97-CD1c_RP_1_1_CD1c`, GSM1717156=`SM_28-CD1c_RP_1_1_CD1c`, 
         GSM1717160=`SM_D36-COK_RP_1_1_CD1c`) %>% 
  ensembl_to_symbol(ensembl) %>% 
  select(-ref_genome) %>% 
  mutate(symbol=transcript, .keep="unused") %>%
  pivot_longer(-c(ensembl, symbol), names_to="sample", values_to="count") %>% 
  nest(data=-sample) %>% 
  unnest(data) %>% 
  mutate(dataset="GSE70106", cell_type="dendritic_myeloid_immature", level=4,
         note="CD11c+ myeloid DC from PBMC")


GSE70106 <- GSE70106_mDC




# 26 Dataset 


# Tally of samples=======================================================================================
cells <- c("dendritic_myeloid_immature", "t_CD4_memory_central", "t_CD4_memory_effector", "dendritic_myeloid_mature",
               "t_CD8_naive", "mast_cell", "t_CD4_memory", "t_helper_h1", "t_helper_h17", "t_helper_h2", "t_reg",
               "t_CD8_memory_central", "macrophage_M2", "t_gamma_delta")
exist_num <- c(2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6, 7)
new_num <- c(20, 33, 39, 0, 13, 13, 24, 14, 15, 25, 35, 15, 5, 4)

data_type <- rep(c("existing", "new"), times=length(cells))
cell_type <- rep(cells, each=2)
count <- seq(2*length(exist_num))
i <- 2*seq_along(exist_num)
count[i] <- new_num
count[i-1] <- exist_num

tibble(cell_type, count, data_type) %>% 
  ggplot(aes(x=cell_type, y=count, fill=data_type, label=count)) +
  geom_col(position = position_stack(reverse = T)) +
  geom_text(size = 3, colour="white", position = position_stack(reverse=T, vjust = 0.5)) +
  labs(y="sample number") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 10, linetype="dashed", colour="blue")

# Combine all data into one dataframe========================================================================

new_data <- bind_rows(GSE135390, GSE138603, GSE138604, GSE122941, GSE113891, GSE85294, GSE89404, GSE125887,
                      GSE115103, GSE151079, GSE107981, GSE75011, GSE118974, GSE155715, GSE147394, GSE136200,
                      GSE118829, GSE123812, GSE153104, GSE133527, GSE131792, GSE115736, GSE70106) %>% 
  mutate(cell_type=as.factor(cell_type))

save(new_data, file = "database_2.RData")


# saveRDS() << correct way
# For the future, save your database as
# 
# 
# 
# Database_table %>%
#   select(sample, symbol, count, dataset, cell_type, note) %>%
#   
#   saveRDS(“counts_second_db_raw.rds”, compression=”gzip”)
# 
# 
# 
# The quotes are to be reformatted

View(new_data)

# GSE135390, GSE113891, GSE107981, GSE75011, GSE118974, GSE136200, GSE118829 don't have ensembl ID
# GSE89404 has NCBI ID instead of ensembl but has symbol. (Solved)
# GSE155715, GSE133527 have ENTREZ ID instead of ENSEMBL ID. (Solved)
# GSE125887 has 3 count columns.
