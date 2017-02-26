# Chapter 5: Ordination in reduced space
# ======================================
# 导入本章所需的程序包 
library(ade4)
library(vegan)
library(gclus)  
library(ape)
# 导入CSV文件数据 
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# 删除没有数据的样方8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]
# 显示环境变量数据集的内容
summary(env)	# 描述性统计
# 全部环境变量数据的PCA分析（基于相关矩阵：参数scale=TURE）
# ********************************************************
env.pca <- rda(env, scale=TRUE) # 参数scale=TRUE 表示对变量进行标准化
env.pca
summary(env.pca) # 默认scaling=2
summary(env.pca, scaling=1)
#注意函数summary（）内的参数scaling，为绘制排序图所选择的标尺类型，
#与函数rda（）内数据标准化的参数scale无关。
# 查看和绘制PCA输出的部分结果 
# ****************************
?cca.object  # 解释vegan包输出的排序结果对象结构和如何提取部分结果
# 特征根
(ev <- env.pca$CA$eig)
# 应用Kaiser-Guttman准则选取排序轴 
ev[ev > mean(ev)]
# 断棍模型
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {
	bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n
bsm
# 绘制每轴的特征根和方差百分比 
par(mfrow=c(2,1))
barplot(ev, main="特征根", col="bisque", las=2)
abline(h=mean(ev), col="red")	# 特征根平均值
legend("topright", "平均特征根", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
	main="% 变差", col=c("bisque",2), las=2)
legend("topright", c("% 特征根", "断棍模型"), 
	pch=15, col=c("bisque",2), bty="n")
	#两种方法是否保留相同的轴数？
# 用一个简单函数绘制与图5.1相同的图：
# 绘制每轴的特征根和方差百分比

source("evplot.R")
evplot(ev)
# 两种PCA双序图：1型标尺和2型标尺 
#********************************
# 使用biplot（）函数绘制排序图 
par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA-1型标尺")
biplot(env.pca, main="PCA-2型标尺")  # 默认 scaling = 2
# 使用cleanplot.pca（）函数绘图
source("cleanplot.pca.R")
cleanplot.pca(env.pca, point=TRUE) 
# point=TRUE表示样方用点表示，变量用箭头表示
cleanplot.pca(env.pca)              
# 默认样方仅用序号标识（同vegan包的标准）
cleanplot.pca(env.pca, ahead=0)    
# ahead=0表示变量用无箭头的线表示
#左图中圆圈代表什么意义？看下面解释内容


# 组合聚类分析结果和排序结果
# ***************************
# 使用环境变量数据对样方进行基于变量标准化后欧氏距离的Ward聚类分析
env.w <- hclust(dist(scale(env)), "ward")
# 裁剪聚类树，只保留4个聚类簇
gr <- cutree(env.w, k=4)
grl <- levels(factor(gr))
# 提取样方坐标，1型标尺
sit.sc1 <- scores(env.pca, display="wa", scaling=1)
# 按照聚类分析的结果对样方进行标识和标色（1型标尺）
p <- plot(env.pca, display="wa", scaling=1, type="n", 
	main="PCA（基于相关矩阵）+聚类簇")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grl)) {
	points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
	}
text(sit.sc1, row.names(env), cex=.7, pos=3)
# 在排序图内添加聚类树
ordicluster(p, env.w, col="dark grey")
legend(locator(1), paste("组",c(1:length(grl))), pch=14+c(1:length(grl)), 
  col=1+c(1:length(grl)), pt.cex=2)
  
# 鱼类多度数据的PCA分析
# **********************
# 物种数据Hellinger转化
spe.h <- decostand(spe, "hellinger")
spe.h.pca <- rda(spe.h)
spe.h.pca
#绘制每轴的特征根和方差百分比
ev <- spe.h.pca$CA$eig
evplot(ev)
# PCA双序图
cleanplot.pca(spe.h.pca, ahead=0)
#这里物种不像环境变量那样能够明显分组。但也能看出物种沿着梯度更替
#分布。在1型标尺的双序图内，可以观察到有8个物种对于第一、二轴有很大
#贡献。可以比较一下，这些物种与4.10.3节聚类分析中聚类簇的指示种是否
#重合？
# 使用PCA（）和biplot.PCA（）两个函数进行环境数据的PCA分析
# **********************************************************
source("PCA.R")  #导入PCA.R脚本，此脚本必须在当前工作目录下或给路径
# PCA，这里默认是1型标尺双序图
env.PCA.PL1 <- PCA(env, stand=TRUE)
biplot.PCA(env.PCA.PL1)
abline(h=0, lty=3)
abline(v=0, lty=3)
# PCA，生成2型标尺双序图
env.PCA.PL2 <- PCA(env, stand=TRUE)
biplot.PCA(env.PCA.PL2 ,scaling=2)
abline(h=0, lty=3)
abline(v=0, lty=3)
#这里主成分轴正负方向是随机的，可能与vegan包输出的排序图成镜像关
#系。但没有关系，因为对象或变量之间的相对位置没有变化。


# 原始物种多度数据的对应分析(CA)
# *******************************
# 计算CA
spe.ca <- cca(spe)
spe.ca
summary(spe.ca)		#默认scaling= 2
summary(spe.ca, scaling=1)
#第一轴有一个很大的特征根。在CA里面，如果特征根超过0.6，代表数据结
#构梯度明显。第一轴特征根占总惯量多少比例呢？需要注意的是，两类标
#尺下，特征根一样。标尺的选择，只影响特征向量，不影响特征根。
# 绘制每轴的特征根和方差百分比
(ev2 <- spe.ca$CA$eig)
evplot(ev2)
#这里，断棍模型比Kaiser-Guttman准则更保守。无论是数量分析结果、还是
#条形图都显示第一轴占绝对优势。
# CA双序图
# *********
par(mfrow=c(1,2))
# 1型标尺：样方点是物种点的形心
plot(spe.ca, scaling=1, main="鱼类多度CA双序图（1型标尺）")
# 2型标尺（默认）：物种点是样方点的形心
plot(spe.ca, main="鱼类多度CA双序图（2型标尺）")
# CA排序中被动加入环境因子
# 调用最后生成CA结果对象（2型标尺）
spe.ca.env <- envfit(spe.ca, env)
plot(spe.ca.env)
# 这个命令的目的是在最后双序图加入环境变量
#新加入的环境变量信息对解读双序图是否有帮助？
# 基于CA排序结果的数据表格重排
# ****************************
vegemite(spe, spe.ca)
#当前输出的表格与传统的群落数据表格排列方式相反，现在以行为物种，以#列为样方。物种排列顺序和样方排列顺序依赖于排序轴的方向（其实是任
#意的）。可以发现，单纯基于第一轴的结果重新排列数据表格，并没有达
#到最佳的效果。因为第二轴所反映的上游（样方1-10）到中游（样方11-18）#的梯度，以及这些样方的特征种，在这个表格里并没有聚集，而是分散的。
# 使用函数CA（）进行对应分析
# ************************
source("CA.R") #导入CA.R脚本，此脚本必须在当前工作目录下或给路径
spe.CA.PL <- CA(spe)
biplot(spe.CA.PL, cex=1)
# 用CA第一轴排序结果重新排列数据表格
# 重新排列数据表格与vegemite（）输出的结果一样
summary(spe.CA.PL)
t(spe[order(spe.CA.PL$F[,1]),order(spe.CA.PL$V[,1])])



# 基于鱼类物种数据Bray-Curtis相异矩阵的PCoA分析
# *********************************************
spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k=(nrow(spe)-1), eig=TRUE)
# 绘制样方主坐标排序图并用加权平均方法将物种投影到样方PCoA排序图
ordiplot(scores(spe.b.pcoa)[,c(1,2)], type="t", main="PCoA分析（带物种投影）") 
abline(h=0, lty=3)
abline(v=0, lty=3)
# 添加物种
spe.wa <- wascores(spe.b.pcoa$points[,1:2], spe)
text(spe.wa, rownames(spe.wa), cex=0.7, col="red")
# 使用pcoa（）函数运行PCoA分析和物种向量投影
# *****************************************
spe.h.pcoa <- pcoa(dist(spe.h))
# 双序图
par(mfrow=c(1,2))
# 第一个双序图：被动加入Hellinger转化的物种数据
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis2=-1) 
abline(h=0, lty=3)
abline(v=0, lty=3)
# 第二个双序图：被动加入Hellinger转化后标准化的物种数据
spe.std <- apply(spe.h, 2, scale)          
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis2=-1) 
abline(h=0, lty=3)
abline(v=0, lty=3)
#如何比较当前PCoA结果与PCA结果？
# 基于欧氏和非欧氏距离的PCoA结果比较
# ***********************************
# 基于Hellinger距离矩阵PCoA
is.euclid(dist(spe.h))
summary(spe.h.pcoa) 
spe.h.pcoa$values
# 基于Bray-Curtis相异矩阵的PCoA
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray) 
spe.bray.pcoa$values    # 观察第18轴及之后的特征根
# 基于Bray-Curtis相异矩阵平方根的PCoA
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values  # 观察特征根
# 基于Bray-Curtis相异矩阵的PCoA（Lingoes校正） 
spe.brayl.pcoa <- pcoa(spe.bray, correction="lingoes")
spe.brayl.pcoa$values   # 观察特征根
# 基于Bray-Curtis相异矩阵的PCoA（Cailliez校正）
spe.brayc.pcoa <- pcoa(spe.bray, correction="cailliez")
spe.brayc.pcoa$values   # 观察特征根
#如果要选择承载最大比例变差的前两轴去了解数据的结构，你会选择上面哪
#种结果呢？



# 基于鱼类数据Bray-Curtis距离矩阵的NMDS排序
# ****************************************
spe.nmds <- metaMDS(spe, distance="bray")
spe.nmds
spe.nmds$stress
plot(spe.nmds, type="t", main=paste("NMDS/Bray-应力函数值=", round(spe.nmds$stress,3)))
#当前生成的排序图与PCA、CA和PCoA的排序图进行比较，有什么不同？
# 评估NMDS拟合度的Shepard图
# **************************
par(mfrow=c(1,2))
stressplot(spe.nmds, main="Shepard图")
gof = goodness(spe.nmds)
plot(spe.nmds, type="t", main="拟合度")
points(spe.nmds, display="sites", cex=gof*200)
# 在NMDS排序图内添加聚类分析结果
# *********************************
# 基于Bray-Curtis相异矩阵Ward聚类结果（提取4组）
spe.bray.ward <- hclust(spe.bray, "ward")
spe.bw.groups <- cutree(spe.bray.ward, k=4)
grp.lev <- levels(factor(spe.bw.groups))
# 与NMDS结果进行组合
sit.sc <- scores(spe.nmds)
p <- ordiplot(sit.sc, type="n", main="NMDS/Bray + clusters Ward/Bray")
for (i in 1:length(grp.lev)) {
	points(sit.sc[spe.bw.groups==i,], pch=(14+i), cex=2, col=i+1)
	}
text(sit.sc, row.names(spe), pos=4, cex=0.7)
# 添加聚类树
ordicluster(p, spe.bray.ward, col="dark grey")
legend(locator(1), paste("Group",c(1:length(grp.lev))), 
  pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=2)
# 一个简单的计算PCA的函数
myPCA <- function(Y) {
Y.mat <- as.matrix(Y)
object.names <- rownames(Y)
var.names <- colnames(Y)
# 中心化数据（计算F矩阵所需）
Y.cent <- scale(Y.mat, center=TRUE, scale=FALSE)
# 协方差矩阵
Y.cov <- cov(Y.cent)
# S的特征向量和特征根（公式9.1和公式9.1）
Y.eig <- eigen(Y.cov)
# 将特征向量赋予矩阵U（用于代表1型标尺双序图变量）
U <- Y.eig$vectors
rownames(U) <- var.names
# 计算F矩阵（用于代表1型标尺双序图对象）
F <- Y.cent%*%U       # 公式9.4
rownames(F) <- object.names
# 计算矩阵U2（用于代表2型标尺双序图变量，见Legendre和Legendre，1998，
# 第397页无标号公式） 
U2 <- U%*%diag(Y.eig$values^0.5)
rownames(U2) <- var.names
# 计算矩阵G（用于代表2型标尺双序图对象，见Legendre和Legendre，1998，第# 404页无标号的公式)
G <- F%*%diag(Y.eig$values^0.5)
rownames(G) <- object.names
# 输出包含所有结果的列表
result <- list(Y.eig$values,U,F,U2,G)
names(result) <- c("eigenvalues","U", "F", "U2", "G")
result
}
# 使用自写函数进行鱼类PCA分析
fish.PCA <- myPCA(spe.h)	
summary(fish.PCA)
# 特征根
fish.PCA$eigenvalues
# 以百分比方式表示特征根
round(100*fish.PCA$eigenvalues/sum(fish.PCA$eigenvalues),2)
# 用总变差（分母）代替特征根和
round(100*fish.PCA$eigenvalues/sum(diag(cov(spe.h))),2)
# 以百分比方式表示累计特征根
round(cumsum(100*fish.PCA$eigenvalues/sum(fish.PCA$eigenvalues)),2)
# 双序图
par(mfrow=c(1,2))
# 1型标尺双序图
biplot(fish.PCA$F, fish.PCA$U)
# 2型标尺双序图
biplot(fish.PCA$G, fish.PCA$U2)


