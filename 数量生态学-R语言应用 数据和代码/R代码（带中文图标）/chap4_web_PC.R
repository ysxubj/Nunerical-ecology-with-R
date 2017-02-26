# ************ Chapter 4 - Cluster analysis *************
# 加载所需的程序包
library(ade4)
library(vegan)  # 应该先加载ade4后加载vegan以避免冲突
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(mvpart)
library(MVPARTwrap) # MVPARTwrap这个程序包必须从本地zip文件安装
# 导入CSV格式的数据
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# 删除无物种数据的样方8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]
# 物种多度数据：先计算样方之间的弦距离矩阵，然后进行单连接聚合聚类
# **************************************************************
spe.norm <- decostand(spe, "normalize")
spe.ch <- vegdist(spe.norm, "euc")
spe.ch.single <- hclust(spe.ch, method="single")
# 使用默认参数选项绘制聚类树
plot(spe.ch.single,main="聚类树",ylab="高度",xlab="单连接聚合聚类")
#基于第一次聚类的结果，如何描述这个数据集？是简单的单一梯度还是区
#分明显的样方组？能否辨认样方的连接链？样方1、5和9为什么最后连
#接？
# 计算完全连接聚合聚类
# ********************
spe.ch.complete <- hclust(spe.ch, method="complete")
plot(spe.ch.complete,main="聚类树",ylab="高度",xlab="完全连接聚合聚类")
#当前所给的样方是沿着河流分布（样方的编号按照流向编排），这个聚类分
#析结果是否将位置相近的样方排在同一个组呢？
#两种完全有效的聚类分析方法分析同一数据，为什么产生如此不同的聚类
#结果呢？
# 计算UPGMA聚合聚类
# ***********************
spe.ch.UPGMA <- hclust(spe.ch, method="average")
plot(spe.ch.UPGMA,main="聚类树",ylab="高度",xlab="UPGMA聚合聚类")
#这个UPGMA聚合聚类树看起来介于单连接聚类和完全连接聚类之间。这种
#情况经常发生。
# 计算鱼类数据的形心聚类
# ***********************
spe.ch.centroid <- hclust(spe.ch, method="centroid")
plot(spe.ch.centroid,main="聚类树",ylab="高度",xlab="UPGMC聚合聚类")
#这种聚类树对生态学家来说简直是噩梦。Legendre 和Legendre（1998，
#第341页）解释了聚类树层级顺序倒置如何产生，并建议用多分法
#（polychotomies）代替二分法（dichotomies）解读这种图。
# 计算Ward最小方差聚类
# ***********************
spe.ch.ward <- hclust(spe.ch, method="ward")
plot(spe.ch.ward,main="聚类树",ylab="高度",xlab="Ward聚类")
#使用距离平方造成此聚类树上半部分过于膨胀。为了使聚类树比例看起来
#更协调而不影响结构，可以使用当前融合水平的平方根重新绘图（图4.5）
spe.ch.ward$height <- sqrt(spe.ch.ward$height)
plot(spe.ch.ward,main="聚类树",ylab="高度",xlab="Ward聚类")

# 同表型相关
# ***********
# 单连接聚类同表型相关
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)
# 完全连接聚类同表型相关
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)
# 平均聚类同表型相关
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)
# Ward聚类同表型相关
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)
#哪个聚类树保持与原始的弦距离矩阵最接近的关系？
#同表型相关也可以用spearman秩相关或Kendall秩相关表示
cor(spe.ch, spe.ch.ward.coph, method="spearman")
# Shepard图
# ***********
par(mfrow=c(2,2))
plot(spe.ch, spe.ch.single.coph, xlab="弦距离", 
  ylab="同表型距离", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)), 
  main=c("单连接",paste("同表型相关",
  round(cor(spe.ch, spe.ch.single.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.single.coph), col="red")
plot(spe.ch, spe.ch.comp.coph, xlab="弦距离", 
	ylab="同表型距离", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
  main=c("完全连接", paste("同表型相关",
  round(cor(spe.ch, spe.ch.comp.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.comp.coph), col="red")
plot(spe.ch, spe.ch.UPGMA.coph, xlab="弦距离", 
	ylab="同表型距离", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)), 
  main=c("UPGMA", paste("同表型相关",
  round(cor(spe.ch, spe.ch.UPGMA.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.UPGMA.coph), col="red")
plot(spe.ch, spe.ch.ward.coph, xlab="弦距离", 
	ylab="同表型距离", asp=1, xlim=c(0,sqrt(2)), 
  ylim=c(0,max(spe.ch.ward$height)),
  main=c("Ward聚类", paste("同表型相关",
  round(cor(spe.ch, spe.ch.ward.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.ward.coph), col="red")
# Gower（1983）距离
gow.dist.single <- sum((spe.ch-spe.ch.single.coph)^2)
gow.dist.comp <- sum((spe.ch-spe.ch.comp.coph)^2)
gow.dist.UPGMA <- sum((spe.ch-spe.ch.UPGMA.coph)^2)
gow.dist.ward <- sum((spe.ch-spe.ch.ward.coph)^2)
gow.dist.single
gow.dist.comp
gow.dist.UPGMA
gow.dist.ward
# 融合水平值图
# *************
par(mfrow=c(2,2))
# 绘制单连接聚类融合水平值图
summary(spe.ch.single)  # 总结聚类分析的结果
plot(spe.ch.single$height, nrow(spe):2, type="S", main="融合水平值-弦距离-单连接", ylab="k （聚类簇数量）", xlab="h (节点高度）)", col="grey")
text(spe.ch.single$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
# 绘制完全连接聚类融合水平值图
plot(spe.ch.complete$height, nrow(spe):2, type="S", 
	main="融合水平值-弦距离-完全连接", 
	ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(spe.ch.complete$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
#在这里建议的组数可能不同。再次裁剪聚类树。这些解决方案是否有意义？
# 绘制UPGMA聚类融合水平值图
plot(spe.ch.UPGMA$height, nrow(spe):2, type="S", 
	main="融合水平值-弦距离-UPGMA", 
	ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(spe.ch.UPGMA$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
#绘制Ward聚类融合水平值图
plot(spe.ch.ward$height, nrow(spe):2, type="S", 
	main="融合水平值-弦距离-Ward", 
	ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(spe.ch.ward$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
#上面4个图看起来差异很大。记住，这些解决方案中任何一个都不是绝对
#正确，每个方案都可以为数据分组提供独特见解。


# 裁剪聚类树以获得k个分类组并使用列联表比较组之间的差异
# *******************************************************
# 设定聚类组的数量
k <- 4  # 根据上面4个融合水平值图，可以观察到分4组水平在所有图里有小
#的跳跃
# 裁剪聚类树
spebc.single.g <- cutree(spe.ch.single, k)
spebc.complete.g <- cutree(spe.ch.complete, k)
spebc.UPGMA.g <- cutree(spe.ch.UPGMA, k)
spebc.ward.g <- cutree(spe.ch.ward, k)
# 通过列联表比较分类结果
# 单连接vs完全连接
table(spebc.single.g, spebc.complete.g)
# 单连接vs UPGMA 
table(spebc.single.g, spebc.UPGMA.g)
#单连接vs Ward 
table(spebc.single.g, spebc.ward.g)
# 完全连接vs UPGMA
table(spebc.complete.g, spebc.UPGMA.g)
# 完全连接vs Ward
table(spebc.complete.g, spebc.ward.g)
# UPGMA vs Ward
table(spebc.UPGMA.g, spebc.ward.g)
#如果两个聚类的结果完全一样，那么这个列联表每行和每列只有一个非零
#数字，其他应该为0。此处并没有出现这种情况。如何解读这些列联表呢？
#例如，单连接聚类第2组含有26个样方，这些样方在Ward聚类中被分散
#到4个组里。


# 依据轮廓宽度图选择最优化的聚类簇数量（Rousseeuw质量指数）
# ********************************************************
# 绘制所有分类水平（除了k=1组的情况）轮廓宽度值（Ward 聚类）
# 首先产生一个长度等于样方数量的空向量asw
asw <- numeric(nrow(spe))
# 其次循环获得轮廓宽度值并依次填入asw向量
for (k in 2:(nrow(spe)-1)) {
	sil <- silhouette(cutree(spe.ch.ward, k=k), spe.ch)
	asw[k] <- summary(sil)$avg.width
	}
# 选择最佳（最大）轮廓宽度值
k.best <- which.max(asw)
# 利用cluster程序包内函数plot.silhouette（）绘制轮廓宽度值k
plot(1:nrow(spe), asw, type="h", 
  main="轮廓宽度-最优聚类簇数（Ward聚类）", 
	xlab="k (组数）", ylab="平均轮廓宽度")
axis(1, k.best, paste("最优",k.best,sep="\n"), col="red", font=2,
  col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
	"with an average silhouette width of", max(asw), "\n")
# 屏幕上将输出如下内容：
# 轮廓宽度最优的聚类簇数量k=2，此时平均轮廓宽度为0.365819
#轮廓宽度法经常选择2组作为最优的分类数量。然而，从生态学角度分析
#分4组似乎更合理。


# 依据Mantel统计（Pearson）相关选择最优分组数量
# **********************************************
# 编写计算代表分类水平的二元距离矩阵的函数
grpdist <- function(X)
  {
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
  }
# 基于Ward 聚类结果运行该函数
kt <- data.frame(k=1:nrow(spe), r=0)
for (i in 2:(nrow(spe)-1)) {
	gr <- cutree(spe.ch.ward, i)
	distgr <- grpdist(gr)
	mt <- cor(spe.ch, distgr, method="pearson")
	kt[i,2] <- mt
}
kt
k.best <- which.max(kt$r)
# 通过cluster程序包内plot.silhouette函数绘制分析图
plot(kt$k, kt$r, type="h", main="Mantel-最优聚类簇数（Ward聚类）", 
	xlab="k (组数）", ylab="Pearson 相关")
axis(1, k.best, paste("最优", k.best, sep="\n"), col="red", font=2,
	col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
# 最终分组的轮廓图
# ****************
# 选择聚类簇的数量
k <- 4
# 轮廓图
cutg <- cutree(spe.ch.ward, k=k)
sil <- silhouette(cutg, spe.ch)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(spe)[attr(silo,"iOrd")]
plot(silo, main="轮廓宽度图-弦距离（Ward聚类）", 
	cex.names=0.8, col=cutg+1, nmax.lab=100,border="white"
)
#组1和组3最连贯，同时组2可能含有被错分的对象。


# 设定组数的最终聚类树
# *********************
# 函数reorder.hclust（）的作用是重新排列从函数hclust（）获得的聚类树，使
# 聚类树内对象的排列顺序与原始相异矩阵内对象的排列顺序尽可能一致。重排# 不影响聚类树的结构。
spe.chwo <- reorder.hclust(spe.ch.ward, spe.ch)
# 绘制重排后带组标识的聚类树
plot(spe.chwo, hang=-1, xlab="4 groups", sub="", ylab="Height", 
	main="Chord - Ward (reordered)", labels=cutree(spe.chwo, k=k))
rect.hclust(spe.chwo, k=k)
# 绘制带不同颜色的最终聚类树 
# 使用我们自编函数hcoplot（）可以快速获得最终聚类树:
source("hcoplot.R")	       # hcoplot.R脚本必须在当前工作目前下
hcoplot(spe.ch.ward, spe.ch, k=4)


# 4个Ward聚类簇在Doub河的分布情况
# ********************************
# 绘制Doubs河流地图（也见第2章）
plot(spa, asp=1, type="n", main="4个Ward聚类组", 
	xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
text(70, 10, "上游", cex=1.2, col="red")
text(15, 115, "下游", cex=1.2, col="red")
# 添加4组分类信息
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for (i in 1:k) {
   points(spa[grw==i,1], spa[grw==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
   }
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("组",1:k), pch=(1:k)+20, col=2:(k+1), 
  pt.bg=2:(k+1), pt.cex=1.5, bty="n")
#请比较当前生成的地图与第2章生成的4种鱼类物种分布地图。



# 热图（heat map）
# ***************
# 用聚类结果重排距离矩阵的热图
dend <- as.dendrogram(spe.chwo)
heatmap(as.matrix(spe.ch), Rowv=dend, symm=TRUE, margin=c(3,3))
#观察如何设定使最热的色彩（计算机输出的是黑色或红色）代表最近的相似
#性，例如对角线代表对象自身的相似性，所以颜色最深。
# 重排群落表格
# 物种按照在样方得分加权平均进行排列
or <- vegemite(spe, spe.chwo)
# 基于聚类树的双排列群落表格的热图
heatmap(t(spe[rev(or$species)]), Rowv=NA, Colv=dend,
	col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
	ylab="物种（样方的加权平均）", xlab="样方")
	
	
	
# 预转化后物种数据k-均值划分
# ****************************
spe.kmeans <- kmeans(spe.norm, centers=4, nstart=100)
#注意：即使给定的nstart相同，每次运行上述命令，所产生的结果也不一定
#完全相同，因为每次运算设定的初始结构是随机的。
# 比较当前分4组的分类结果与之前Ward聚类的结果。
table(spe.kmeans$cluster, spebc.ward.g) 
#这两个聚类结果是否非常相似？哪个（或哪些）对象有差别？


# k-均值划分，2组到10组
# ************************
spe.KM.cascade <- cascadeKM(spe.norm, inf.gr=2, sup.gr=10, iter=100, 
  criterion="ssi")
plot(spe.KM.cascade, sortg=TRUE,border="white")
#该图显示每个对象在每种分类组数下的归属（图上每行代表一种组数）。图
#内的表格有不同的颜色，每行两种颜色，代表分两组k=2，三种颜色代表k=3，#依此类推。右图代表不同k值条件下的中止标准的统计量。此系列中，到底
#多少组数是最佳方案？如果倾向于较大的组数，哪个是最佳方案呢？

summary(spe.KM.cascade)
spe.KM.cascade$results
#这里，最小SSE值用于确定在给定组数k下的最佳方案，而calinski和ssi
#指标用于确定最佳k值，两个指标解决不同的问题。
spe.KM.cascade$partition
#记住不同的组数k={2,3,…，10}，每次都是独立运行。因此图4.15从下到上
#结构互相独立，与层次聚类树的嵌套结构不同。
# 按照k-均值划分结果重新排列样方
spe[order(spe.kmeans$cluster),]
# 使用函数vegemite（）重新排列样方-物种矩阵
ord.KM <- vegemite(spe, spe.kmeans$cluster)
spe[ord.KM$sites, ord.KM$species]


# 形心点划分（PAM）
# 基于弦距离矩阵进行计算
# ***********************
# 聚类簇数量的选择
asw <- numeric(nrow(spe))
# 循环计算2至28组分类数下平均轮廓宽度
# asw对象名称取自"average silhouette width"
for (k in 2:(nrow(spe)-1)) 
	asw[k] <- pam(spe.ch, k, diss=TRUE)$silinfo$avg.width
k.best <- which.max(asw) 
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
	"with an average silhouette width of", max(asw), "\n")
plot(1:nrow(spe), asw, type="h", main="Choice of the number of clusters", 
	xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
	col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
#当k=2时，PAM具有最佳的方案（asw=0.3841），这并不是我们期望的结果。
#如果选择常用的4组分法，从轮廓宽度角度分析结果并不好（asw=0.2736）。
#尽管如此，我们还是需要分析PAM分4组的情况。
# PAM分4组情况
spe.ch.pam <- pam(spe.ch, k=4, diss=TRUE)
summary(spe.ch.pam)
spe.ch.pam.g <- spe.ch.pam$clustering
spe.ch.pam$silinfo$widths
# 将当前的分类结果与之前Ward聚类和k-均值划分进行比较
table(spe.ch.pam.g, spebc.ward.g)
table(spe.ch.pam.g, spe.kmeans$cluster)
#PAM结果与Ward聚类和k-均值划分结果显著不同
# k=4组下k-均值法和PAM法轮廓宽带图
par(mfrow=c(1,2))
plot(silhouette(spe.kmeans$cluster,spe.ch), main="轮廓宽度图：k-均值法", 
  cex.names=0.8, col= spe.kmeans$cluster+1,border="white")
plot(silhouette(spe.ch.pam), main="轮廓宽度图：PAM", cex.names=0.8, 
	col=spe.ch.pam$silinfo$widths+1,border="white")
#基于此图，可以分辨PAM和k-均值法哪个有更好的平均轮廓值。请将当前
#结果与Ward聚类比较轮廓宽带图。


# 绘制4组k-均值样方分组结果沿Doubs河的分布地图
# ***********************************************
plot(spa, asp=1, type="n", main="Four k-means groups", 
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
text(50, 10, "Upstream", cex=1.2, col="red")
text(25, 115, "Downstream", cex=1.2, col="red")
grKM <- spe.kmeans$cluster
k <- length(levels(factor(grKM)))
for (i in 1:k) {
   points(spa[grKM==i,1], spa[grKM==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
   }
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("Group",1:k), pch=(1:k)+20, col=2:(k+1), 
  pt.bg=2:(k+1), pt.cex=2, bty="n")
  
  
# 基于k-均值划分结果（分4组）的样方聚类簇与4个环境变量的关系
# *************************************************************
attach(env)
# 定量环境变量箱线图
# 海拔、坡度、氧含量和铵浓度
par(mfrow=c(2,2))
boxplot(sqrt(alt) ~ spe.kmeans$cluster, main="海拔", las=1, 
  ylab="sqrt(alt)", col=2:5, varwidth=TRUE)
boxplot(log(pen) ~ spe.kmeans$cluster, main="坡度", las=1, 
  ylab="log(pen)", col=2:5, varwidth=TRUE)
boxplot(oxy ~ spe.kmeans$cluster, main="氧含量", las=1, 
  ylab="oxy", col=2:5, varwidth=TRUE)
boxplot(sqrt(amm) ~ spe.kmeans$cluster, main="铵浓度", las=1, 
  ylab="sqrt(amm)", col=2:5, varwidth=TRUE)
# 检验方差分析的假设
# 检验残差的正态性
shapiro.test(resid(lm(sqrt(alt) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(log(pen) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(oxy ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(sqrt(amm) ~ as.factor(spe.kmeans$cluster))))
#检验结果表明sqrt(alt)、log(pen)、oxy和sqrt(amm)的残差是正态分布。
#也尝试为其他的环境变量寻找好的标准化转化。
# 检验方差齐性
bartlett.test(sqrt(alt), as.factor(spe.kmeans$cluster))
bartlett.test(log(pen), as.factor(spe.kmeans$cluster))
bartlett.test(oxy, as.factor(spe.kmeans$cluster))
bartlett.test(sqrt(amm), as.factor(spe.kmeans$cluster))
#变量sqrt(alt)的方差不齐，所以参数检验的方差分析不适用
#可被参数检验的变量的方差分析 
summary(aov(log(pen) ~ as.factor(spe.kmeans$cluster)))
summary(aov(oxy ~ as.factor(spe.kmeans$cluster)))
summary(aov(sqrt(amm) ~ as.factor(spe.kmeans$cluster)))
#坡度、含氧量和铵浓度在不同聚类簇之间是否差异显著？
# 变量alt的Kruskal-Wallis检验
kruskal.test(alt ~ as.factor(spe.kmeans$cluster))
#海拔在不同聚类簇之间是否有显著差异？
detach(env)


 # 两种聚类的列联表 
# ******************
# 基于环境变量（见第2章）的样方聚类
env2 <- env[,-1]
env.de <- vegdist(scale(env2), "euc")
env.kmeans <- kmeans(env.de, centers=4, nstart=100)
env.KM.4 <- env.kmeans$cluster
# 比较两种聚类的结果
table(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))
#两种聚类结果是否相同？
# 用卡方检验分析两种聚类之间的差异 
chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))
# 改用置换的方法进行卡方检验分析 
chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster), 
	simulate.p.value=TRUE)
	
	
# k-均值样方聚类簇平均多度
# ***********************
groups <- as.factor(spe.kmeans$cluster)
spe.means <- matrix(0, ncol(spe), length(levels(groups)))
row.names(spe.means) <- colnames(spe)
for(i in 1:ncol(spe)) {
  spe.means[i,] <- tapply(spe[,i], spe.kmeans$cluster, mean)
  }
# 每组的物种平均多度
group1 <- round(sort(spe.means[,1], decreasing=TRUE), 2)
group2 <- round(sort(spe.means[,2], decreasing=TRUE), 2)
group3 <- round(sort(spe.means[,3], decreasing=TRUE), 2)
group4 <- round(sort(spe.means[,4], decreasing=TRUE), 2)
# 显示多度大于平均值的物种
group1.domin <- which(group1>mean(group1))
group1
group1.domin 

#...对其他组进行相同的分析
# Kendall一致性系数（W）
# ********************
# 提取多度最大的物种
sp.sum <- apply(spe, 2, sum)
spe.sorted <- spe[,order(sp.sum, decreasing=TRUE)]
spe.small <- spe.sorted[,1:20]
# 物种数据转化和矩阵转置
spe.small.hel <- decostand(spe.small, "hellinger")
spe.small.std <- decostand(spe.small.hel, "standardize")
spe.small.t <- t(spe.small.std)
# 物种k-均值划分
spe.t.kmeans.casc <- cascadeKM(spe.small.t, inf.gr=2, sup.gr=8,
	iter=100, criterion="calinski")
plot(spe.t.kmeans.casc, sortg=TRUE)
#结果显示分2组是不错的选择，如果分更多的组，可能有些组只有一个物种。
# 在object $partition中显示分为2组
clusters <- spe.t.kmeans.casc$partition[,1]
clusters
#现在有两个物种组，对这两组运行全局的Kendall W检验
spe.kendall.global <- kendall.global(spe.small.hel, clusters)
spe.kendall.global
spe.kendall.post <- kendall.post(spe.small.hel, clusters, nperm=9999)
spe.kendall.post

# 基于有-无数据的物种集合
# ***********************
source("test.a.R")                # 函数必须在当前工作目录下
spe.pa <- decostand(spe,"pa")      # 将数据转化为有-无数据
res <- test.a(spe.pa, nperm=99999)
summary(res)
#输出的结果是res$p.a.dsit包含一个p值的矩阵。下一步计算向量化的p值矩阵
#的Holm校正数（见7.2.6节）。
res.p.vec <- as.vector(res$p.a.dist)
adjust.res <- p.adjust(res.p.vec,method="holm")
#在Holm校正后的p值中，可以发现0.05或比0.05小但很接近的值。例如，数
#值0.04878的未校正p值应该是0.00018。Holm校正后的p值矩阵内大约有83
#个值小于0.05（可能每次运行这个数值有轻微变化）。
#因此0.00018或更小概率才能称为显著。接下来将相异矩阵内将所有大于#0.00018的值都换成1。
res.pa.dist <- res$p.a.dist
res.pa.dist[res.pa.dist>0.00018] <- 1
#由p值构成的相异矩阵可以用热图表示（见第3章）
source("coldiss.R")               # 函数必须在当前工作目录下
coldiss(res.pa.dist, nc=10, byrank=TRUE, diag=TRUE)
# 运行模糊非层次聚类（见4.12节）。尝试使用多k值获得最佳的分组。
res.pa.fuz <- fanny(res.pa.dist, k=5, memb.exp=1.5)
summary(res.pa.fuz)
plot(silhouette(res.pa.fuz), main="Silhouette plot - Fuzzy clustering", cex.names=0.8, col=res.pa.fuz$silinfo$widths+1)
res.pa.fuz$silinfo                 # 轮廓信息
#此轮廓信息包含一个聚类簇中最可能的成员（"cluster"）、最近邻体聚类
#簇（"neighbor"）和轮廓宽度值("sil_width")。





# Dufrene & Legendre物种指示值
# ****************************
# 依据das（离源头距离）环境变量将样方分为4组
das.D1 <- dist(data.frame(das=env[,1], row.names=rownames(env)))
dasD1.kmeans <- kmeans(das.D1, centers=4, nstart=100)
dasD1.kmeans$cluster
# 样方组的指示种
(iva <- indval(spe, dasD1.kmeans$cluster) )
#输出结果包括以下表格：
#-relfrg=物种在每个组的相对频度
#-relabu=物种在组间的相对多度
#-indval=每个物种的指示值
#下面两项内容将从indval表格内提取出来（含有最高指示值的组）和置换检
#验的结果。
# 显著指示物种的表格
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(spe>0, 2, sum)[iva$pval<=0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg
# 将结果输出为CSV文件（可以用电子表格打开）
write.csv(fidg, "IndVal-das.csv")
#基于上面的结果，你能获得什么信息？这个结果是否与第1章所描述的4种
#生态区域有关联？


# 多元回归树
# **********
library(mvpart) # 载入mvpart程序包
spe.ch.mvpart <- mvpart(data.matrix(spe.norm) ~ ., env, margin=0.08, cp=0, 
	xv="pick", xval=nrow(spe), xvmult=100, which=4)
# 此处可以点击想要的分组数（例如4组）
summary(spe.ch.mvpart)
printcp(spe.ch.mvpart)
# MRT的残差
par(mfrow=c(1,2))
hist(residuals(spe.ch.mvpart), col="grey")
plot(predict(spe.ch.mvpart), residuals(spe.ch.mvpart), 
	main="Residuals vs Predicted")
abline(h=0, lty=3, col="grey")
# 组的组成
spe.ch.mvpart$where
# 识别组的名称
(groups.mrt <- levels(as.factor(spe.ch.mvpart$where)))
# 第一片叶子的鱼类物种组成 
spe.norm[which(spe.ch.mvpart$where==groups.mrt[1]),]
# 第一片叶子的环境变量组成
env[which(spe.ch.mvpart$where==groups.mrt[1]),]
# 叶子的鱼类物种组成表格和饼图
leaf.sum <- matrix(0, length(groups.mrt), ncol(spe))
colnames(leaf.sum) <- colnames(spe)
for(i in 1:length(groups.mrt)){
  leaf.sum[i,] <- apply(spe.norm[which(spe.ch.mvpart$where==groups.mrt[i]),],
    2, sum)
}
leaf.sum
par(mfrow=c(2,2))
for(i in 1:length(groups.mrt)){
	pie(which(leaf.sum[i,]>0), radius=1, main=c("leaf #", groups.mrt[i]))
}
# 从mvpart（）函数获得的结果对象中提取MRT结果
# 必须加载MVPARTwrap和rdaTest程序包
spe.ch.mvpart.wrap <- MRT(spe.ch.mvpart, percent=10, species=colnames(spe))
summary(spe.ch.mvpart.wrap)


# 在MRT的结果中寻找指示种
spe.ch.MRT.indval <- indval(spe.norm, spe.ch.mvpart$where)
spe.ch.MRT.indval$pval    # 概率
# 为每个显著的物种寻找最高指示值的叶子spe.ch.MRT.indval$maxcls[which(spe.ch.MRT.indval$pval<=0.05)]
# 每个显著的物种在最高指示值的叶子中的指示值
spe.ch.MRT.indval$indcls[which(spe.ch.MRT.indval$pval<=0.05)]
# MRT作为空间和时间系列约束聚类方法 
spe.ch.seq <- mvpart(as.matrix(spe) ~ das, env, cp=0, xv="pick", margin=0.08,
  xval=nrow(spe), xvmult=100, which=4)
# 此时可以点击所期望的分组组数
summary(spe.ch.seq)
# 组的组成（终端节点的标识）
(gr <- spe.ch.seq$where)
# 重新编排聚类簇的编号
aa <- 1
gr2 <- rep(1,length(gr))
for (i in 2:length(gr)) {
  if (gr[i]!=gr[i-1]) aa <- aa+1
  gr2[i] <- aa
}
# 在Doubs河地图上标注样方所属的聚类簇
plot(spa, asp=1, type="n", main="MRT 分类组", 
	xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
text(70, 10, "上游", cex=1.2, col="red")
text(15, 115, "下游", cex=1.2, col="red")
k <- length(levels(factor(gr2)))
for (i in 1:k) {
   points(spa[gr2==i,1], spa[gr2==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
   }
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("组",1:k), pch=(1:k)+20, col=2:(k+1), 
  pt.bg=2:(k+1), pt.cex=2, bty="n")


# 鱼类数据的c-均值模糊聚类
# ************************
k <- 4  # 选择聚类分组的数量
spe.fuz <- fanny(spe.ch, k=k, memb.exp=1.5)
summary(spe.fuz)
spefuz.g <- spe.fuz$clustering
# 样方成员值
spe.fuz$membership
# 每个样方最接近的聚类簇
spe.fuz$clustering
# 轮廓图
plot(silhouette(spe.fuz), main="鱼类数据的c-均值模糊聚类", 
  cex.names=0.8, col=spe.fuz$silinfo$widths+1,border="white")
# 模糊聚类簇的主坐标排序(PCoA)
dc.pcoa <- cmdscale(spe.ch)
dc.scores <- scores(dc.pcoa, choices=c(1,2))
plot(scores(dc.pcoa), asp=1, type="n",
	main="模糊聚类簇的排序(PCoA)")
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
for (i in 1:k) {
	gg <- dc.scores[spefuz.g==i,]
	hpts <- chull(gg)
	hpts <- c(hpts, hpts[1])
	lines(gg[hpts,], col=i+1)
	}
stars(spe.fuz$membership, location=scores(dc.pcoa), draw.segments=TRUE,
	add=TRUE, scale=FALSE, len=0.1, col.segments=2:(k+1))
legend(locator(1), paste("组", 1:k, sep=" "),
	pch=15, pt.cex=2, col=2:(k+1), bty="n")
# 在图上任意点击某一处放置图例
