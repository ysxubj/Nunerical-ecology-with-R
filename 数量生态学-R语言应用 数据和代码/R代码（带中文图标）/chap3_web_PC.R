# 加载所需程序包
library(ade4)
library(vegan)	# 应该先加载ade4再加载vegan，以避免一些冲突
library(gclus)
library(cluster)
library(FD)
# 导入CSV格式的数据
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# 剔除无物种数据的样方8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]
# 定量（半定量）数据的相异和距离测度
# **********************************
# 原始物种数据的Bray-Curtis相异矩阵
spe.db <- vegdist(spe)	# Bray-Curtis相异系数（默认）
head(spe.db)
# 对数转化后物种数据的Bray-Curtis相异矩阵
spe.dbln <- vegdist(log1p(spe))
head(spe.dbln)
# 弦距离矩阵
spe.norm <- decostand(spe, "nor")
spe.dc <- dist(spe.norm)
head(spe.dc)
# Hellinger距离矩阵
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
head(spe.dh)
#查看decostand（）函数和vegdist（）函数的帮助文件，寻找计算卡方距离#矩阵的参数设定，卡方距离对当前这些数据也适用。
# 二元数据的相异测度
# ******************
#注意：所有的二元距离函数在计算系数时，均会自动对数据进行二元转化
#因此这里的数据不需要二元转化（decostand(，"pa")）。函数dist.binary（）
#会自动对数据进行二元转化，但函数vegist（）需要设定参数binary=TRUE。
# 使用vegdist（）函数计算Jaccard相异矩阵
spe.dj <- vegdist(spe, "jac", binary=TRUE)
head(spe.dj)
head(sqrt(spe.dj))
# 使用dist（）函数计算Jaccard相异矩阵
spe.dj2 <- dist(spe, "binary")
head(spe.dj2)
# 使用dist.binary（）函数计算Jaccard相异矩阵
spe.dj3 <- dist.binary(spe, method=1)
head(spe.dj3)
# 使用dist.binary（）函数计算S?rensen相异矩阵
spe.ds <- dist.binary(spe, method=5)
head(spe.ds)
# 使用vegdist（）函数计算S?rensen相异矩阵
spe.ds2 <- vegdist(spe, binary=TRUE)
head(spe.ds2)
head(sqrt(spe.ds2))
# Ochiai相异矩阵
spe.och <- dist.binary(spe, method=7)
head(spe.och)
#这里显示Jaccard和S?rensen距离矩阵有两种数值。可以回看3.3节引言部
#分了解这两种数值的差异。
# 图解关联矩阵
# ************
# 运行coldiss（）函数时会自动要求加载gclus程序包，我们也可以预先加载
library(gclus)
# 导入我们自己编写的函数coldiss（）
source("coldiss.R")  #如果函数文件没有在当前工作目录下，需要指定文件路径
# 使用coldiss（）函数产生的彩图（在一些数据分析文献中也称作热图或格状图）
# ******************************************************************
#绘制两种（原始和重排）相异矩阵图并标色。屏幕输出颜色分类：红紫色=
#相异系数接近0（最大相似系数）；青绿色=相异系数接近1（最小相似系数）
#coldiss（）函数使用说明：
# coldiss(D=dist.object, nc=4, byrank=TRUE, diag=FALSE)
# D应该是一个相异矩阵
#如果D为最大值大于1的距离矩阵，此时D会除以max(D)
# nc颜色种类数量
# byrank=TRUE   等大小分级，即每个颜色所包含的值的数量一样多
# byrank=FALSE  等区间分级，即每个颜色所包含的值的区间一样长
#如果diag=TRUE，表示样方号放置在矩阵对角线上
# 图解基于原始鱼类多度数据的Bray-Curtis相异矩阵
# 等区间分级的4种颜色（方便比较）
coldiss(spe.db, byrank=FALSE, diag=TRUE)
# 图解基于对数转化数据的Bray-Curtis相异矩阵
coldiss(spe.dbln, byrank=FALSE, diag=TRUE)
# 弦距离矩阵
coldiss(spe.dc, byrank=FALSE, diag=TRUE)
# Hellinger距离矩阵
coldiss(spe.dh, byrank=FALSE, diag=TRUE)
# Jaccard距离矩阵
coldiss(spe.dj, byrank=FALSE, diag=TRUE)
#请比较当前的Jaccard距离热图和前面的其他距离矩阵热图。Jaccard图是
#基于二元数据计算。这是否影响结果呢？Jaccard图和前面的数量系数图之
#间的差异是否比数量系数图之间的差异大呢？
# 简单匹配相异系数（在ade4程序包内也称为Sokal & Michener指数)
spe.s1 <- dist.binary(spe, method=2)
coldiss(spe.s1^2, byrank=FALSE, diag=TRUE)
#比较一下当前的对称相异矩阵与之前的Jaccard矩阵。哪个受双零问题影响
#更突出？
# 剔除env数据框内das变量
env2 <- env[,-1]
# 由标准化后的env2数据框计算的欧氏距离矩阵
env.de <- dist(scale(env2))
coldiss(env.de, diag=TRUE)
# 物种数据的Hellinger距离矩阵（等数量的分级）
coldiss(spe.dh, diag=TRUE)
#因为具有相同的样方位置排序，可以比较图3.2的左边热图与前面基于环
#境变量的热图。你是否能观察到一些共同的特征？
# 基于二维空间坐标的欧氏距离矩阵
spa.de <- dist(spa)
coldiss(spa.de, diag=TRUE)
# 基于一维das变量（离源头距离）的欧氏距离矩阵
das.df <- as.data.frame(env$das, row.names=rownames(env))
riv.de <- dist(das.df)
coldiss(riv.de, diag=TRUE)
#为什么基于x-y的欧氏距离图和基于das的欧氏距离图有这样的差异？
# 生成30个对象、5个二元变量的数据集，每个变量有预先设置的固定的0和1的
# 数量
# 变量1：10个1和20个0，顺序随机
var1 <- sample(c(rep(1,10), rep(0,20)))
# 变量2：15个0在一起，15个1在一起
var2 <- c(rep(0,15), rep(1,15))
# 变量3：3个1，3个0交替出现，直到总数量达到30为止
var3 <- rep(c(1,1,1,0,0,0),5)
# 变量4：5个1，10个0交替出现，直到总数量达到30为止
var4 <- rep(c(rep(1,5), rep(0,10)), 2)
# 变量5：前16元素是7个1和9个0的随机排列，接着是4个0和10个1
var5.1 <- sample(c(rep(1,7), rep(0,9)))
var5.2 <- c(rep(0,4), rep(1,10))
var5 <- c(var5.1, var5.2)
# 将变量1至变量5合成一个数据框
dat <- data.frame(var1, var2, var3, var4, var5)
dim(dat)
# 简单匹配系数的计算（在ade4程序包中也称为Sokal & Michener指数）
dat.s1 <- dist.binary(dat, method=2)
coldiss(dat.s1, diag=TRUE)
# 为计算Gower指数（S15）的虚拟数据
# 随机生成30个平均值为0、标准差为1的正态分布数据
var.g1 <- rnorm(30,0,1)
# 随机生成30个从0到5均匀分布的数据
var.g2 <- runif(30,0,5)
# 生成3个水平的因子变量（每个水平10个重复）
var.g3 <- gl(3,10)
# 生成与var.g3正交的2个水平的因子变量
var.g4 <- gl(2,5,30)
#var.g3和var.g4组合在一起代表一个双因素交叉平衡设计
dat2 <- data.frame(var.g1,var.g2,var.g3,var.g4)
summary(dat2)
# 使用daisy（）函数计算Gower相异矩阵
# 完整的Gower相异矩阵（基于4个变量）
dat2.S15 <- daisy(dat2, "gower")
range(dat2.S15)
coldiss(dat2.S15, diag=TRUE)
# 仅使用2个正交的因子变量计算gower相异矩阵
dat2partial.S15 <- daisy(dat2[,3:4], "gower")
coldiss(dat2partial.S15, diag=TRUE)
# 在dat2partial.S15矩阵内相异系数值代表什么？
levels(factor(dat2partial.S15))
#对象对所对应的数值分享共同的因子水平2、1和无因子水平。最高相异系数
#值的对象对不分享共同的因子水平。
# 使用FD程序包内gowdis（）函数计算Gower相异矩阵
library(FD)  #如果FD还没载入
?gowdis
dat2.S15.2 <- gowdis(dat2)
range(dat2.S15.2)
coldiss(dat2.S15.2, diag=TRUE)
# 仅使用两个正交的因子变量计算距离矩阵
dat2partial.S15.2 <- gowdis(dat2[,3:4])
coldiss(dat2partial.S15.2, diag=TRUE)
# 在dat2partial.S15.2矩阵内相异系数值代表什么？
levels(factor(dat2partial.S15.2))

# R模式相异矩阵
# *************
# 物种多度矩阵的转置矩阵
spe.t <- t(spe)
# 先卡方转化后计算欧氏距离
spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)
coldiss(spe.t.D16, diag=TRUE)
   #在右边的图中，你能否分辨出物种组？
# 鱼类有-无数据的Jaccard指数
spe.t.S7 <- vegdist(spe.t, "jaccard", binary=TRUE)
coldiss(spe.t.S7, diag=TRUE)
#将右边的图与之前获得的卡方距离图进行比较，物种组是否一致？

# 环境变量之间的Pearson线性相关系数r
env.pearson <- cor(env)	# 默认 method = "pearson"
round(env.pearson, 2)
# 在绘图之前重新排位变量
env.o <- order.single(env.pearson)
# pairs：一个同时生成双变量之间散点图和相关系数图的函数
# 图中上半部分显示两个变量之间相关系数（带显著水平)
source("panelutils.R") # 如果脚本不在当前工作目录下，需要给出访问路径
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
	diag.panel=panel.hist, main="Pearson 相关矩阵")
par(op)
#请辨认与变量"das"相关的环境变量。这些图能够告诉你什么信息呢？
# 环境变量之间的Kendall秩相关
env.ken <- cor(env, method="kendall")
env.o <- order.single(env.ken)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
	method="kendall", diag.panel=panel.hist, main="Kendall Correlation Matrix")
par(op)
#通过这些双变量关系图，你更倾向于用Kendall秩相关还是Pearson线性相关？
