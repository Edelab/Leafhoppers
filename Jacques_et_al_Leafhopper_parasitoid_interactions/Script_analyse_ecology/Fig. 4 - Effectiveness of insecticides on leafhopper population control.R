# Pesticides

# Pesticide with at least 4 applications was selected

# Libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggeffects)
library(esquisse)
library(lubridate)
library(readxl)
library(ggpubr)

# Load the data
g1<-read_excel("pest.xlsx", 2)
g1$AI<-as.factor(g1$AI)
g1$Species<-as.factor(g1$Species)
g1$Time<-factor(g1$Time, levels=c("Before", "After"))

# Subset per species
Efa<-subset(g1, Species=="Empoasca fabae")
Mqa<-subset(g1, Species=="Macrosteles quadrilineatus")
Locsp<-subset(g1, Species=="Local species")


# Base graphic
esquisser()

# Empoasca----

efatest_results <- Efa %>%
  group_by(AI) %>%
  do({
    before <- filter(., Time == "Before")$Number
    after <- filter(., Time == "After")$Number
    test_output <- wilcox.test(before, after, paired = TRUE)
    data.frame(
      p.value = test_output$p.value,
      V.statistic = test_output$statistic
    )
  })

# View the results
print(efatest_results)

efa_results<-as.data.frame(efatest_results)

write.xlsx(efa_results, file = "efa_results.xlsx")


Efa_fig <- ggplot(Efa) +
  aes(x = Time, y = Number, fill = Time) +
  geom_boxplot() +
  scale_fill_manual(values = c(Before = "#7fbf5b", After="#d6eac7")) +
  labs(x = "Time", y = "Number of insects normalized") +
  theme_classic(12) +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(vars(AI)) +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.signif", 
                     label.y = max(Efa$Number) * 0.8,
                     label.x = 1.5,
                     size=6,
                     color="black")

Efa_fig


ggsave("Efa_fig.tiff", width = 180, height = 180, units = "mm",dpi=600, compression = "lzw")
ggsave("Efa_fig.png", width = 180, height = 180, units = "mm",dpi=600)



# Macrosteles----

mqatest_results <- Mqa %>%
  group_by(AI) %>%
  do({
    before <- filter(., Time == "Before")$Number
    after <- filter(., Time == "After")$Number
    test_output <- wilcox.test(before, after, paired = TRUE)
    data.frame(
      p.value = test_output$p.value,
      V.statistic = test_output$statistic
    )
  })

# View the results
print(mqatest_results)

mqa_results<-as.data.frame(mqatest_results)

write.xlsx(mqa_results, file = "mqa_results.xlsx")


Mqa_fig <- ggplot(Mqa) +
  aes(x = Time, y = Number, fill = Time) +
  geom_boxplot() +
  scale_fill_manual(values = c(Before = "#e27a1c", After="#ffd2b3")) +
  labs(x = "Time", y = "Number of insects normalized") +
  theme_classic(12) +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(vars(AI)) +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.signif", 
                     label.y = max(Mqa$Number) * 0.8,
                     label.x = 1.5,
                     size=6,
                     color="black")

Mqa_fig


ggsave("Mqa_fig.tiff", width = 180, height = 180, units = "mm",dpi=600, compression = "lzw")
ggsave("Mqa_fig.png", width = 180, height = 180, units = "mm",dpi=600)


# Local----

loctest_results <- Locsp %>%
  group_by(AI) %>%
  do({
    before <- filter(., Time == "Before")$Number
    after <- filter(., Time == "After")$Number
    test_output <- wilcox.test(before, after, paired = TRUE)
    data.frame(
      p.value = test_output$p.value,
      V.statistic = test_output$statistic
    )
  })

# View the results
print(loctest_results)

loctest_results<-as.data.frame(loctest_results)

write.xlsx(loctest_results, file = "loctest_results.xlsx")


Loc_fig <- ggplot(Locsp) +
  aes(x = Time, y = Number, fill = Time) +
  geom_boxplot() +
  scale_fill_manual(values = c(Before = "#1f74b4", After="#dfe6f3")) +
  labs(x = "Time", y = "Number of insects normalized") +
  theme_classic(12) +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(vars(AI)) +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.signif", 
                     label.y = 15,
                     label.x = 1.5,
                     size=6,
                     color="black")

Loc_fig


ggsave("Loc_fig.tiff", width = 180, height = 180, units = "mm",dpi=600, compression = "lzw")
ggsave("Loc_fig.png", width = 180, height = 180, units = "mm",dpi=600)



# Paper figure
g1<-read_excel("pest.xlsx", 4)
g1$AI<-as.factor(g1$AI)
g1$Species<-factor(g1$Species,levels=c("Empoasca fabae",
                                       "Macrosteles quadrilineatus",
                                       "Local species"))
g1$Time<-factor(g1$Time, levels=c("Before", "After"))
g1$Timecol<-factor(g1$Timecol, levels=c("Efab", "Efaa",
                                        "Mqab", "Mqaa",
                                        "Lspb", "Lspa"))



Figpaper <- ggplot(g1) +
  aes(x = Time, y = Number) +
  geom_boxplot(aes(fill=Timecol)) +
  scale_fill_manual(values = c( Efab = "#7fbf5b", Efaa="#d6eac7",
                                Mqab = "#e27a1c", Mqaa="#ffd2b3", 
                               Lspb = "#1f74b4", Lspa="#dfe6f3")) +
  labs(x = "Time", y = "Number of insects normalized") +
  theme_classic(12) +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  facet_grid(Species~AI, scales="free_y") +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.signif", 
                     label.y = 15,
                     label.x = 1.5,
                     size=5,
                     color="black")

Figpaper

ggsave("Figpaper.tiff", width = 160, height = 185, units = "mm",dpi=600, compression = "lzw")
ggsave("Figpaper.png", width = 160, height = 185, units = "mm",dpi=600)
