library(tidyverse)
library(ggplot2)

# Load data
df1 <- read.csv("data/snps.metrics.tsv", header = TRUE, na.strings = c("."), sep ="\t")
df2 <- read.csv("data/indels.metrics.tsv", header = TRUE, na.strings = c("."), sep ="\t")



# Replace the dots by NA


# Organize data
df1 <- df1 %>% mutate(type = 'indels')
df2 <- df2 %>% mutate(type = 'snps')
df3 <- bind_rows(df1,df2)
df3 <- df3 %>% mutate_if(is.factor, as.numeric)
df3 <- df3 %>% mutate_if(is.integer, as.numeric)

# Plots annotation...

# 1 - FS density
ggplot(data=df3) + geom_density(aes(x=FS, group=type, fill=type), alpha=0.5) +
  scale_x_continuous(limits=c(0, 20))
p <- ggplot(data=df3) + geom_density(aes(x=FS, group=type, fill=type), alpha=0.5) +
                 scale_x_continuous(limits=c(0, 20))
ggsave(plot=p, paste("FS.svg", sep=""), width = 14, height = 7, dpi=600)

# 2 - SOR density
p <- ggplot(data=df3) + geom_density(aes(x=SOR, group=type, fill=type), alpha=0.5) +
                 scale_x_continuous(limits=c(0, 5))
ggsave(plot=p, paste("SOR.svg", sep=""), width = 14, height = 7, dpi=600)



# 3 - MQRankSum density
p <- ggplot(data=df3) + geom_density(aes(x=MQRankSum, group=type, fill=type), alpha=0.5) +
                  scale_x_continuous(limits=c(-2, 4))
ggsave(plot=p, paste("MQRanksum.svg", sep=""), width = 14, height = 7, dpi=600)


# 4 - ReadPosRankSum
p <- ggplot(data=df3) + geom_density(aes(x=ReadPosRankSum, group=type, fill=type), alpha=0.5) +
                  scale_x_continuous(limits=c(-5, 20))
ggsave(plot=p, paste("ReadPosRankSum.svg", sep=""), width = 14, height = 7, dpi=600)

# 5 - QD density
p <- ggplot(data=df3) + geom_density(aes(x=QD, group=type, fill=type), alpha=0.5)
ggsave(plot=p, paste("QD.svg", sep=""), width = 14, height = 7, dpi=600)

# 6 - MQ density
p <- ggplot(data=df3) + geom_density(aes(x=MQ, group=type, fill=type), alpha=0.5,adjust=0.001) +
                  scale_x_continuous(limits=c(59, 61))
ggsave(plot=p, paste("MQ.svg", sep=""), width = 14, height = 7, dpi=600)

# 7 - DP density
p <- ggplot(data=df3) + geom_density(aes(x=DP, group=type, fill=type), alpha=0.5) +
                 scale_x_continuous(limits=c(0, 200))
ggsave(plot=p, paste("DP.svg", sep=""), width = 14, height = 7, dpi=600)


# 8 - QUAL density
p <- ggplot(data=df3) + geom_density(aes(x=QUAL, group=type, fill=type), alpha=0.5) + 
                  scale_x_continuous(limits=c(0, 200))
ggsave(plot=p, paste("QUAL.svg", sep=""), width = 14, height = 7, dpi=600)



