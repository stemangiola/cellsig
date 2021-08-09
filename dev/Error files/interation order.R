# This file shows that the order of interaction terms may affect the summary plot

saveRDS(final_tSNE, "final_tSNE.rds")

final_tSNE <- readRDS("dev/intermediate_data/final_tSNE.rds")

analysis_sec <- final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5)) +
  ggtitle("interaction(ancestor_type, analysis)")

ggsave("analysis_sec.png", analysis_sec)

analysis_fir <- final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(analysis, ancestor_type), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5)) +
  ggtitle("interaction(analysis, ancestor_type)")

ggsave("analysis_fir.png", analysis_fir)
