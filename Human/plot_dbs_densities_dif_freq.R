library(ggplot2)
library(grid)
library(gridExtra)

args=commandArgs(TRUE)

pathogenic_filename=args[1]
likely_path_filename=args[2]
benign_filename=args[3]
ref_pathogenic_filename=args[4]
ref_likely_path_filename=args[5]
ref_benign_filename=args[6]
benign_freq_filename=args[7]

pathogenic_bss <- read.table(pathogenic_filename, col.names=c("sequence_name", "score"))
pathogenic_bss$type <- "pathogenic"
pathogenic_bss$freq <- "pathogenic/\nlikely pathogenic"
likely_path_bss <- read.table(likely_path_filename, col.names=c("sequence_name", "score"))
likely_path_bss$type <- "pathogenic"
likely_path_bss$freq <- "pathogenic/\nlikely pathogenic"
benign_bss <- read.table(benign_filename, col.names=c("sequence_name", "score"))
benign_bss$type <- "benign"
benign_bss_freq <- read.table(benign_freq_filename, col.names=c("sequence_name", "freq"))
benign_bss <- merge(benign_bss, benign_bss_freq, by="sequence_name")
benign_bss$freq <- as.character(benign_bss$freq)
benign_bss$freq[benign_bss$freq=="10"]<-"benign 10%"
benign_bss$freq[benign_bss$freq=="20"]<-"benign 20%"
benign_bss$freq[benign_bss$freq=="50"]<-"benign 50%"

all_bss <- rbind(pathogenic_bss, likely_path_bss, benign_bss)

ref_pathogenic_bss <- read.table(ref_pathogenic_filename, col.names=c("sequence_name", "refscore"))
ref_likely_path_bss <- read.table(ref_likely_path_filename, col.names=c("sequence_name", "refscore"))
ref_benign_bss <- read.table(ref_benign_filename, col.names=c("sequence_name", "refscore"))
all_ref_bss <-rbind(ref_pathogenic_bss, ref_likely_path_bss, ref_benign_bss)

all_dbss <- merge(all_bss,all_ref_bss,by="sequence_name")

all_dbss$DBscore <- all_dbss$refscore-all_dbss$score

all_dbss_freq <- all_dbss[!(all_dbss$type=="benign" & all_dbss$freq == "clinVar"), ]

pdf("plots_dif_freq/DBscores.pdf")

ggplot(all_dbss_freq, aes(x=freq, y=DBscore, colour=freq))+geom_point()

dev.off()

pdf("plots_dif_freq/DBscore_dencities.pdf")

ggplot(all_dbss_freq, aes(x=DBscore, fill=freq)) + geom_density(alpha=0.25)+ scale_fill_manual(values=c("blueviolet", "blue", "green", "red"))

dev.off()

all_dbss_freq$freq[all_dbss_freq$freq=="benign 10%"]<-"benign 10-20%"
all_dbss_freq$freq[all_dbss_freq$freq=="benign 20%"]<-"benign 10-20%"

p1 <- ggplot(all_dbss_freq,aes(x=DBscore, fill=freq)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("olivedrab3","deepskyblue3", "red3"))+scale_size_manual(values=c(10,5))+ theme(legend.position="none")+labs(x="", y="")+theme(axis.text=element_text(size=12))
p2 <- ggplot(all_dbss_freq,aes(x=DBscore, fill=freq)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("olivedrab3","deepskyblue3", "red3"),name = "Frequencies\n")+xlim(-5, 10)+labs(y="", x="Delat-bitscore")+ theme(legend.text=element_text(size=10), legend.title=element_text(size=14))+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

plots <- list(p1, p2)

tleft <- textGrob("Density", gp=gpar(fontsize = 16), rot = 90)
ttop <- textGrob("DBS of pathogenic/likely pathogenic and benign\nvariants with different minor allele frequencies\n", gp=gpar(fontsize = 22))

pdf("plots_dif_freq/DBscore_dencities_10-20_combined.pdf")

grid.arrange(grobs=plots, ncol=1, left=tleft, top=ttop)

dev.off()

#dbs_path=all_dbss$DBscore[all_dbss$type=="pathogenic"]
#dbs_ben=all_dbss$DBscore[all_dbss$type=="benign"]
#
#pdf("plots_likely-path/DBscore_dencities.pdf")
#
#plot(density(dbs_ben), col="blue", main="DBscore", xlim=c(min(c(dbs_ben, dbs_path)), max(c(dbs_ben, dbs_path))) )
#lines(density(dbs_path), col="red")
#legend(20, 0.18, legend=c("pathogenic/likely pathogenic", "benign"),
#       col=c("red", "blue"), lty=1, cex=1.3)
#
#dev.off()


#mylogit <- glm(type ~ DBscore, data=all_dbss, family=binomial())






