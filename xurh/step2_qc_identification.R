# step2_quality control on identification result
rm(list = ls())
options(stringsAsFactors = F)
source("udf.R")
creDir("data")
creDir("res")
creDir("public")
dataDir = file.path(getwd(),'data')
resDir = file.path(getwd(),'res')
publicDir = file.path(getwd(),'public')
library(ggplot2)
library(ggThemeAssist)
load("output_step1.RData")
# score statistics
dat <- dat[order(dat$"Score", decreasing = TRUE),]
score_table <- as.data.frame(table(cut(
  dat$Score,
  breaks = c(
    min(dat$Score),
    5,
    10,
    50,
    100,
    500,
    1000,
    2000,
    5000,
    max(dat$Score)
  ),
  include.lowest = TRUE
)))
p1 <- ggplot(data = score_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of protein Sequest.HT score") +
  ylab("Protein count") +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.05*max(score_table$Freq))) +
  scale_x_discrete(expand = c(0, 0))
# ggThemeAssistGadget(p1)
p1 <- p1 + theme(
    axis.line = element_line(linetype = "solid"), 
    axis.text = element_text(face = "bold"), 
    axis.text.x = element_text(hjust = 0.6, vjust = 0.5, 
        angle = 45), panel.background = element_rect(fill = NA))
p1
ggsave(file.path(resDir, "score_bar.jpeg"))

# score distribution plot
p2 <-
  ggplot(data = dat,
         mapping = aes(x = 1:length(dat$"Score"), y = Score)) +
  geom_point() +
  xlab("A point is a protein") +
  ylab("Sequest.HT Score of protein")
p2 <- p2 + theme(
  axis.line = element_line(linetype = "solid"), 
  axis.text = element_text(face = "bold"), 
  panel.background = element_rect(fill = NA))
p2
ggsave(file.path(resDir, "score_dist.jpeg"), dpi = 300 )

# Coverage statistics
summary(dat$Coverage)
coverage_table <- as.data.frame(table(cut(
  dat$Coverage,
  breaks = seq(0, 100, 10),
  include.lowest = TRUE
)))
p3 <- ggplot(data = coverage_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of protein coverage (%)") +
  ylab("Protein count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.05*max(coverage_table$Freq)))
p3 <- p3 + theme(
  axis.line = element_line(linetype = "solid"), 
  axis.text = element_text(face = "bold"), 
  panel.background = element_rect(fill = NA))
p3
ggsave(file.path(resDir, "coverage_bar.jpeg"), dpi = 300 )

# unique peptides statistics
summary(dat$"#.Unique.Peptides")
unipep_table <- as.data.frame(table(cut(
  dat$"#.Unique.Peptides",
  breaks = c(
    min(dat$"#.Unique.Peptides"),
    1,
    5,
    10,
    20,
    30,
    40,
    50,
    100,
    max(dat$"#.Unique.Peptides")
  ),
  include.lowest = TRUE
)))
p4 <- ggplot(data = unipep_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of unique peptides number") +
  ylab("Protein count") +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.05*max(unipep_table$Freq)))
p4 <- p4 + theme(
  axis.line = element_line(linetype = "solid"), 
  axis.text = element_text(face = "bold"), 
  axis.text.x = element_text(hjust = 0.6, vjust = 0.5, 
                             angle = 45), panel.background = element_rect(fill = NA))
p4
ggsave(file.path(resDir, "unipep_bar.jpeg"))

# unique peptides distribution plot
dat <- dat[order(dat$"#.Unique.Peptides" , decreasing = TRUE),]
p5 <- ggplot(data = dat,
         mapping = aes(x = 1:length(dat$"#.Unique.Peptides"), y = dat$"#.Unique.Peptides")) +
  geom_point() +
  xlab("A point is a protein") +
  ylab("Number of unique peptides")
p5 <- p5 + theme(
  axis.line = element_line(linetype = "solid"), 
  axis.text = element_text(face = "bold"), 
  axis.text.x = element_text(hjust = 0.6, vjust = 0.5, 
                             angle = 45), panel.background = element_rect(fill = NA))
p5
ggsave(file.path(resDir, "unipep_dist.jpeg"))
