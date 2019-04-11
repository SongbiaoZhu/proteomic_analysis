# step2_quality control on identification result
rm(list = ls())

source("functions.R")
cre_folder("out_folder")

library(ggplot2)

load("output_step1.RData")

theme_set(theme_gray()+theme(axis.line = element_line(size=0.5),panel.background = element_rect(fill=NA,size=rel(20)), panel.grid.minor = element_line(colour = NA), axis.text = element_text(size=16), axis.title = element_text(size=18)))
# quality control on score
# score statistics
dat <- dat[order(dat$"Score.Sequest.HT", decreasing = TRUE),]
summary(dat$Score.Sequest.HT)
score_table <- as.data.frame(table(cut(
  dat$Score.Sequest.HT,
  breaks = seq(min(dat$Score.Sequest.HT), max(dat$Score.Sequest.HT), 1000),
  include.lowest = TRUE
)))
score_table <- as.data.frame(table(cut(
  dat$Score.Sequest.HT,
  breaks = c(
    min(dat$Score.Sequest.HT),
    5,
    10,
    50,
    100,
    500,
    1000,
    2000,
    5000,
    max(dat$Score.Sequest.HT)
  ),
  include.lowest = TRUE
)))
ggplot(data = score_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of protein Sequest.HT score") +
  ylab("Protein count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1.05*max(score_table$Freq)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5)
ggsave("out_folder/score_bar.jpeg")
# score distribution plot
score_dist <-
  ggplot(data = dat,
         mapping = aes(x = 1:length(dat$"Score.Sequest.HT"), y = Score.Sequest.HT)) +
  geom_point() +
  xlab("A point is a protein") +
  ylab("Sequest.HT Score of protein")
score_dist
ggsave("out_folder/score_dist.jpeg")

# quality control on Coverage
# Coverage statistics
summary(dat$Coverage)
coverage_table <- as.data.frame(table(cut(
  dat$Coverage,
  breaks = seq(min(dat$Coverage), max(dat$Coverage), 10),
  include.lowest = TRUE
)))
ggplot(data = coverage_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of protein coverage (%)") +
  ylab("Protein count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5)
ggsave("out_folder/coverage_bar.jpeg")

# unique peptides statistics
summary(dat$"X..Unique.Peptides")
unipep_table <- as.data.frame(table(cut(
  dat$X..Unique.Peptides,
  breaks = c(
    min(dat$X..Unique.Peptides),
    1,
    5,
    10,
    20,
    30,
    40,
    50,
    100,
    max(dat$X..Unique.Peptides)
  ),
  include.lowest = TRUE
)))
ggplot(data = unipep_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Intervals of unique peptides number") +
  ylab("Protein count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5)
ggsave("out_folder/unipep_bar.jpeg")
# unique peptides distribution plot
dat <- dat[order(dat$X..Unique.Peptides , decreasing = TRUE),]
unipep_dist <-
  ggplot(data = dat,
         mapping = aes(x = 1:length(dat$"X..Unique.Peptides"), y = X..Unique.Peptides)) +
  geom_point() +
  xlab("A point is a protein") +
  ylab("Number of unique peptides")
unipep_dist
ggsave("out_folder/unipep_dist.jpeg")
