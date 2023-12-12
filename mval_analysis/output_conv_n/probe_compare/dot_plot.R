library(tidyverse)

mval <- read.csv("mval_sum2.csv")
ggplot(data=mval) +
  geom_point(mapping=aes(x=norm_avg, y=c1_avg, shape = conv_set, color = conv_set)) + 
  geom_abline(slope = 1 , intercept = 0.4519978172) +
  geom_abline(slope = 1 , intercept = -0.4519978172) +
  geom_abline(slope = 1 , intercept = 0.941069716, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.941069716, linetype = "dashed") +
  labs(title= "cl1 conv. probes deviation from normal tissue m-val",
     caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
# probe must have at least on cluster with a pval < 0.05 to be included
# according to the spreadsheet, using the average change as cut off will reduce the probes by ~38%
ggplot(data=mval) +
  geom_point(mapping=aes(x=norm_avg, y=c2_avg, shape = conv_set, color = conv_set)) + 
  geom_abline(slope = 1 , intercept = 0.4519978172) +
  geom_abline(slope = 1 , intercept = -0.4519978172) +
  geom_abline(slope = 1 , intercept = 0.941069716, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.941069716, linetype = "dashed") +
  labs(title= "cl2 conv. probes deviation from normal tissue m-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
ggplot(data=mval) +
  geom_point(mapping=aes(x=norm_avg, y=c3_avg, shape = conv_set, color = conv_set)) + 
  geom_abline(slope = 1 , intercept = 0.4519978172) +
  geom_abline(slope = 1 , intercept = -0.4519978172) +
  geom_abline(slope = 1 , intercept = 0.941069716, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.941069716, linetype = "dashed") +
  labs(title= "cl3 conv. probes deviation from normal tissue m-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")


bval <- read.csv("bval_sum2.csv")
ggplot(data=bval) +
  geom_point(mapping=aes(x=norm_avg, y=c1_avg, shape = conv_set, color = conv_set)) +
  geom_abline(slope = 1 , intercept = 0.07062241962) +
  geom_abline(slope = 1 , intercept = -0.07062241962) +
  geom_abline(slope = 1 , intercept = 0.1584036328, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.1584036328, linetype = "dashed") +
  labs(title= "cl1 conv. probes deviation from normal tissue B-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
ggplot(data=bval) +
  geom_point(mapping=aes(x=norm_avg, y=c2_avg, shape = conv_set, color = conv_set)) +
  geom_abline(slope = 1 , intercept = 0.07062241962) +
  geom_abline(slope = 1 , intercept = -0.07062241962) +
  geom_abline(slope = 1 , intercept = 0.1584036328, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.1584036328, linetype = "dashed") +
  labs(title= "cl2 conv. probes deviation from normal tissue B-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
ggplot(data=bval) +
  geom_point(mapping=aes(x=norm_avg, y=c3_avg, shape = conv_set, color = conv_set)) +
  geom_abline(slope = 1 , intercept = 0.07062241962) +
  geom_abline(slope = 1 , intercept = -0.07062241962) +
  geom_abline(slope = 1 , intercept = 0.1584036328, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -0.1584036328, linetype = "dashed") +
  labs(title= "cl3 conv. probes deviation from normal tissue B-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) or cl1 (dashed)")

############################################################################

mval_vol <- read.csv("mval_sum.csv")
ggplot(data=mval_vol) +
  geom_point(mapping=aes(x=c1_deviation, y=-log10(pval_norm_vs_c1), shape = conv_set, color = conv_set)) + 
  geom_vline(xintercept=c(0.4519978172,-0.4519978172)) +
  geom_vline(xintercept=c(0.941069716,-0.941069716), linetype = "dashed") +
  geom_hline(yintercept=0.05) +
  labs(title= "Conv. probes deviation from normal tissue M-val",
       caption = "y-int set to p = 0.05; x-int set to the std.dev. of the change in cl3 (line) or cl1 (dashed)")

############################################################################

mval_opt <- read.csv("mval_opt2.csv")
ggplot(data=mval_opt) +
  geom_point(mapping=aes(x=norm_avg, y=c1_avg)) +
  geom_abline(slope = 1 , intercept = 0.4138559439) +
  geom_abline(slope = 1 , intercept = -0.4138559439) +
  geom_abline(slope = 1 , intercept = 1.113764978, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -1.113764978, linetype = "dashed") +
  labs(title= "cl1 conv. probes deviation from normal tissue m-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
# probe must have at least on cluster with a pval < 0.1 to be included
ggplot(data=mval_opt) +
  geom_point(mapping=aes(x=norm_avg, y=c2_avg)) +
  geom_abline(slope = 1 , intercept = 0.4138559439) +
  geom_abline(slope = 1 , intercept = -0.4138559439) +
  geom_abline(slope = 1 , intercept = 1.113764978, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -1.113764978, linetype = "dashed") +
  labs(title= "cl1 conv. probes deviation from normal tissue m-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
# avg change in cl2 deviation is 0.1182836599
ggplot(data=mval_opt) +
  geom_point(mapping=aes(x=norm_avg, y=c3_avg)) +
  geom_abline(slope = 1 , intercept = 0.4138559439) +
  geom_abline(slope = 1 , intercept = -0.4138559439) +
  geom_abline(slope = 1 , intercept = 1.113764978, linetype = "dashed") +
  geom_abline(slope = 1 , intercept = -1.113764978, linetype = "dashed") +
  labs(title= "cl1 conv. probes deviation from normal tissue m-val",
       caption = "slopes set on std.dev. of the change in cl3 (line) of cl1 (dashed)")
# avg change in cl3 deviation is -0.06581890

############################################################################

mval_opt_vol <- read.csv("mval_opt.csv")
ggplot(data=mval_opt_vol) +
  geom_point(mapping=aes(x=c1_deviation, y=-log10(pval_norm_vs_c1))) + 
  xlim(-2.2, 2.2) +
  ylim(-0.1, 13) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1)

ggplot(data=mval_opt_vol) +
  geom_point(mapping=aes(x=c2_deviation, y=-log10(pval_norm_vs_c2))) + 
  xlim(-2.2, 2.2) +
  ylim(-0.1, 13) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1)

ggplot(data=mval_opt_vol) +
  geom_point(mapping=aes(x=c3_deviation, y=-log10(pval_norm_vs_c3))) + 
  xlim(-2.2, 2.2) +
  ylim(-0.1, 13) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1)

############################################################################

mval_opt_sum_vol <- read.csv("mval_opt_sum.csv")
mval_opt_sum_vol$c1_deviation <- mval_opt_sum_vol$c1_avg - mval_opt_sum_vol$norm_avg
mval_opt_sum_vol$c2_deviation <- mval_opt_sum_vol$c2_avg - mval_opt_sum_vol$norm_avg
mval_opt_sum_vol$c3_deviation <- mval_opt_sum_vol$c3_avg - mval_opt_sum_vol$norm_avg
used <- read.csv("used_opt.csv")
mval_opt_sum_vol$probes <- used$e0

ggplot(mval_opt_sum_vol, aes(x=c1_deviation, y=-log10(pval_norm_vs_c1) ) ) +
  geom_hex(bins = 200) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c1_deviation, y=-log10(pval_norm_vs_c1), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1) + 
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.00360")

ggplot(mval_opt_sum_vol, aes(x=c2_deviation, y=-log10(pval_norm_vs_c2) ) ) +
  geom_hex(bins = 200) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c2_deviation, y=-log10(pval_norm_vs_c2), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1) +
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.03310")

ggplot(mval_opt_sum_vol, aes(x=c3_deviation, y=-log10(pval_norm_vs_c3) ) ) +
  geom_hex(bins = 200) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c3_deviation, y=-log10(pval_norm_vs_c3), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439)) +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed") +
  geom_hline(yintercept=0.1) + 
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.22598")

############################################################################

ggplot(mval_opt_sum_vol, aes(x=c1_deviation, y=-log10(pval_norm_vs_c1) ) ) +
  geom_point(data = mval_opt_sum_vol, mapping=aes(x=c1_deviation, y=-log10(pval_norm_vs_c1))) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c1_deviation, y=-log10(pval_norm_vs_c1), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439), color = "blue") +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed", color = "blue") +
  geom_hline(yintercept=0.05) + 
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.00360")

ggplot(mval_opt_sum_vol, aes(x=c2_deviation, y=-log10(pval_norm_vs_c2) ) ) +
  geom_point(data = mval_opt_sum_vol, mapping=aes(x=c2_deviation, y=-log10(pval_norm_vs_c2))) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c2_deviation, y=-log10(pval_norm_vs_c2), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439), color = "blue") +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed", color = "blue") +
  geom_hline(yintercept=0.05) +
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.03310")

ggplot(mval_opt_sum_vol, aes(x=c3_deviation, y=-log10(pval_norm_vs_c3) ) ) +
  geom_point(data = mval_opt_sum_vol, mapping=aes(x=c3_deviation, y=-log10(pval_norm_vs_c3))) +
  geom_point(data = mval_opt_sum_vol[which(mval_opt_sum_vol$probes>0),], mapping=aes(x=c3_deviation, y=-log10(pval_norm_vs_c3), size = 2, color = "pink")) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlim(-3, 3) +
  ylim(-0.1, 15) +
  geom_vline(xintercept=c(0.4138559439,-0.4138559439), color = "blue") +
  geom_vline(xintercept=c(1.113764978,-1.113764978), linetype = "dashed", color = "blue") +
  geom_hline(yintercept=0.05) + 
  labs(caption = "section are divided according to std.dev. of the change in cl3 (line) or cl1 (dashed)
       fisher p-value for significance of the association is 0.22598")

############################################################################


ggplot(data=em_genes, aes(x = fct_reorder(term_description, percent), y=percent)) + 
  ylim(0, 25) +
  labs(title= "Enrichment",
       caption = "percent = number of associated genes on list found in David / all genes listed found in David") + 
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.title = element_text(size = 10), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.5, "lines")) +
  coord_flip()


