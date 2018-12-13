library(ggplot2)
library(grid)
library(gridExtra)
library(pROC)
args=commandArgs(TRUE)

pathogenic_filename=args[1]
likely_path_filename=args[2]
benign_filename=args[3]
ref_pathogenic_filename=args[4]
ref_likely_path_filename=args[5]
ref_benign_filename=args[6]

pathogenic_bss <- read.table(pathogenic_filename, col.names=c("sequence_name", "score"))
pathogenic_bss$type <- "pathogenic"
likely_path_bss <- read.table(likely_path_filename, col.names=c("sequence_name", "score"))
likely_path_bss$type <- "likely pathogenic"
benign_bss <- read.table(benign_filename, col.names=c("sequence_name", "score"))
benign_bss$type <- "benign"

all_bss <- rbind(pathogenic_bss, likely_path_bss, benign_bss)

pdf("plots_likely-path/bitscores.pdf")

ggplot(all_bss, aes(x=type, y=score, colour=type))+geom_point()

dev.off()

ref_pathogenic_bss <- read.table(ref_pathogenic_filename, col.names=c("sequence_name", "refscore"))
ref_likely_path_bss <- read.table(ref_likely_path_filename, col.names=c("sequence_name", "refscore"))
ref_benign_bss <- read.table(ref_benign_filename, col.names=c("sequence_name", "refscore"))
all_ref_bss <-rbind(ref_pathogenic_bss, ref_likely_path_bss, ref_benign_bss)

all_dbss <- merge(all_bss,all_ref_bss,by="sequence_name")

all_dbss$DBscore <- all_dbss$refscore-all_dbss$score

pdf("plots_likely-path/DBscores.pdf")

ggplot(all_dbss, aes(x=type, y=DBscore, colour=type))+geom_point()

dev.off()

all_dbss_only_path <- all_dbss[all_dbss$type=="pathogenic" | all_dbss$type=="benign",]

all_dbss_path_lp_combained <- all_dbss


all_dbss_path_lp_combained$type[all_dbss_path_lp_combained$type=="likely pathogenic"]<-"pathogenic/\nlikely pathogenic"
all_dbss_path_lp_combained$type[all_dbss_path_lp_combained$type=="pathogenic"]<-"pathogenic/\nlikely pathogenic"

pdf("plots_likely-path/DBscore_densities_ggplot.pdf")

ggplot(all_dbss_path_lp_combained,aes(x=DBscore, fill=type)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("deepskyblue3", "red3"),name = "Variants\n")+labs(x="Delat bit-score")+ ggtitle("DBS for pathogenic, likely pathogenic and benign variants")+ coord_fixed(ratio = 100)+theme(plot.title = element_text(hjust = 0.5))

dev.off()

p1 <- ggplot(all_dbss_only_path,aes(x=DBscore, fill=type)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("deepskyblue3", "red3"),name = "Variants\n")+labs(x="", y="")+ coord_fixed(ratio = 100)+theme(axis.text=element_text(size=10))+ theme(legend.text=element_text(size=10), legend.title=element_text(size=14))
p2 <- ggplot(all_dbss_path_lp_combained,aes(x=DBscore, fill=type)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("deepskyblue3","red3"),name = "Variants\n")+labs(x="", y="")+ coord_fixed(ratio = 100)+theme(axis.text=element_text(size=10))+ theme(legend.text=element_text(size=10), legend.title=element_text(size=14))
p3 <- ggplot(all_dbss,aes(x=DBscore, fill=type)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("deepskyblue3", "orange", "red3"),name = "Variants\n")+labs(y="", x="Delta-bitscore")+ coord_fixed(ratio = 100)+theme(axis.text=element_text(size=10), axis.title=element_text(size=16))+ theme(legend.text=element_text(size=10), legend.title=element_text(size=14))

#p1 <- ggplot(all_dbss,aes(x=DBscore, fill=type)) + geom_density(alpha=0.25)+ scale_fill_manual(values=c("blue", "red"))+ theme(legend.position="none")+labs(x="", y="")
#p2 <- ggplot(all_dbss,aes(x=DBscore, fill=type)) + geom_density(alpha=0.25)+ scale_fill_manual(values=c("blue", "red"),name = "Variants\n")+xlim(-5, 10)+labs(y="", x="Delat bit-score")

tleft <- textGrob("Density", gp=gpar(fontsize = 16), rot = 90)
ttop <- textGrob("DBS of pathogenic, likely-pathogenic and benign variants\n", gp=gpar(fontsize = 18))


plots <- list(p1, p2, p3)

pdf("plots_likely-path/DBscore_densities_path_vs_likely_path.pdf")

grid.arrange(grobs=plots, ncol=1, left=tleft, top=ttop)

dev.off()

plot.roc(all_dbss_path_lp_combained$type, all_dbss_path_lp_combained$DBscore, print.thres=TRUE, print.thres.col="red", print.auc=TRUE)



