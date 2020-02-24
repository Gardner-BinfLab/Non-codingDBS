library(ggplot2)

args=commandArgs(TRUE)

pathogenic_filename=args[1]
benign_filename=args[2]
ref_pathogenic_filename=args[3]
ref_benign_filename=args[4]

pathogenic_bss <- read.table(pathogenic_filename, col.names=c("sequence_name", "score"))
pathogenic_bss$type <- "pathogenic"
benign_bss <- read.table(benign_filename, col.names=c("sequence_name", "score"))
benign_bss$type <- "benign"

all_bss <- rbind(pathogenic_bss, benign_bss)

pdf("plots/bitscores.pdf")

ggplot(all_bss, aes(x=type, y=score, colour=type))+geom_point()

dev.off()

ref_pathogenic_bss <- read.table(ref_pathogenic_filename, col.names=c("sequence_name", "refscore"))
ref_benign_bss <- read.table(ref_benign_filename, col.names=c("sequence_name", "refscore"))
all_ref_bss <-rbind(ref_pathogenic_bss,ref_benign_bss)

all_dbss <- merge(all_bss,all_ref_bss,by="sequence_name")

all_dbss$DBscore <- all_dbss$refscore-all_dbss$score

pdf("plots/DBscores.pdf")

ggplot(all_dbss, aes(x=type, y=DBscore, colour=type))+geom_point()

dev.off()


pdf("plots/DBscore_densities_ggplot.pdf")

ggplot(all_dbss,aes(x=DBscore, fill=type)) + geom_density(alpha=0.5)+ scale_fill_manual(values=c("deepskyblue3", "red3"),name = "Variants\n")+labs(x="Delat bit-score")+ ggtitle("DBS for pathogenic and benign variants")+ coord_fixed(ratio = 100)+theme(plot.title = element_text(hjust = 0.5))

dev.off()


dbs_path=all_dbss$DBscore[all_dbss$type=="pathogenic"]
dbs_ben=all_dbss$DBscore[all_dbss$type=="benign"]

pdf("plots/DBscore_dencities.pdf")

plot(density(dbs_ben), col="blue", main="DBscore", xlim=c(min(c(dbs_ben, dbs_path)), max(c(dbs_ben, dbs_path))) )
lines(density(dbs_path), col="red")
legend(20, 0.18, legend=c("pathogenic", "benign"),
       col=c("red", "blue"), lty=1, cex=1.3)

dev.off()


#mylogit <- glm(type ~ DBscore, data=all_dbss)






