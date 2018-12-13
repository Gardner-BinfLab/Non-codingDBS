library(pROC)
library(boot)
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

ref_pathogenic_bss <- read.table(ref_pathogenic_filename, col.names=c("sequence_name", "refscore"))
ref_likely_path_bss <- read.table(ref_likely_path_filename, col.names=c("sequence_name", "refscore"))
ref_benign_bss <- read.table(ref_benign_filename, col.names=c("sequence_name", "refscore"))
all_ref_bss <-rbind(ref_pathogenic_bss, ref_likely_path_bss, ref_benign_bss)

all_dbss <- merge(all_bss,all_ref_bss,by="sequence_name")

all_dbss$DBscore <- all_dbss$refscore-all_dbss$score

all_dbss_only_path <- all_dbss[all_dbss$type=="pathogenic" | all_dbss$type=="benign",]

all_dbss_path_lp_combained <- all_dbss

all_dbss_path_lp_combained$type[all_dbss_path_lp_combained$type=="likely pathogenic"]<-"pathogenic/\nlikely pathogenic"
all_dbss_path_lp_combained$type[all_dbss_path_lp_combained$type=="pathogenic"]<-"pathogenic/\nlikely pathogenic"

all_dbss_path_lp_combained$type <- as.factor(all_dbss_path_lp_combained$type)

mylogit <- glm(type ~ DBscore, data=all_dbss_path_lp_combained, family=binomial())

summary(mylogit)

cv.err=cv.glm(all_dbss_path_lp_combained, mylogit)
cv.err$delta
cv.error=cv.glm(all_dbss_path_lp_combained, mylogit, K=10)$delta
cv.error

predict(mylogit, type="response")

pdf("ROCcurve.pdf")
plot.roc(all_dbss_path_lp_combained$type, all_dbss_path_lp_combained$DBscore, print.thres=TRUE, print.thres.col="red", print.auc=TRUE)
dev.off()


