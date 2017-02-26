# CHAPTER 7: SPATIAL ANALYSIS
# ***************************
# 载入本章所用的程序包
library(ape)
library(spdep)
library(vegan)
library(ade4)
# 以下几个程序包可以从https://r-forge.r-project.org/R/?group_id=195下载
library(packfor)	
library(spacemakeR)	
library(AEM)	
library(PCNM)	
source("plot.links.R")  # 脚本函数必须在当前工作文件夹中
source("sr.value.R")   # 脚本函数必须在当前工作文件夹中
# 导入数据
mite <- read.table("mite.txt")
mite.env <- read.table("mite_env.txt")
mite.xy <- read.table("mite_xy.txt")
mite.h <- decostand (mite, "hellinger")
mite.xy.c <- scale(mite.xy, center=TRUE, scale=FALSE)
# 空间相关图（基于Moran指数)
# ***********************
# 寻找距离在0.7m范围内的所有样方对，并计算样方对之间的滞后阶数
plot.links(mite.xy, thresh=0.7)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, 0.7)
summary(nb1)
# 基质密度的相关图
subs.dens <- mite.env[,1]
subs.correlog <- sp.correlogram(nb1, subs.dens, order=14, method="I", zero.policy=TRUE)
print(subs.correlog, p.adj.method="holm")
plot(subs.correlog)
# 甲螨数据Mantel相关图
# ********************
# 数据首先进行去趋势处理（见第7.3.2节）
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy)) 
mite.h.D1 <- dist(mite.h.det)
(mite.correlog <- mantel.correlog(mite.h.D1, XY=mite.xy, nperm=99))
summary(mite.correlog)   
plot(mite.correlog)
# 等级数
mite.correlog$n.class #或：mite.correlog[2]
# 分割点
mite.correlog$break.pts #或：mite.correlog[3]
# 趋势面分析
# ***********
# 从规则的正方形表面取样的简单模型
# 构建一个10×10的栅格样区
xygrid <- expand.grid(1:10, 1:10)
plot(xygrid)
xygrid.c <- scale(xygrid, scale=FALSE)  # 中心化
X <- xygrid.c[,1]
Y <- xygrid.c[,2]
# 绘制X和Y的一阶、二阶和三阶函数
par(mfrow=c(3,3))
s.value(xygrid,(X))
s.value(xygrid,(Y))
s.value(xygrid,(X+Y))
s.value(xygrid,(X^2+Y^2))
s.value(xygrid,(X^2-X*Y-Y^2))
s.value(xygrid,(X+Y+X^2+X*Y+Y^2))
s.value(xygrid,(X^3+Y^3))
s.value(xygrid,(X^3+X^2*Y+X*Y^2+Y^3))
s.value(xygrid,(X+Y+X^2+X*Y+Y^2+X^3+X^2*Y+X*Y^2+Y^3))
	#可以尝试其他组合，比如使用减号或回归系数不等于1的情况。
# 甲螨数据趋势面分析
# *******************
# 构建中心化甲螨数据X-Y坐标标准（非正交）的三阶多项式函数
mite.poly <- poly(as.matrix(mite.xy.c), degree=3, raw=TRUE)
colnames(mite.poly) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")
#poly（）函数产生多项式按顺序分别是：X，X^2，X^3，Y，XY， X^2Y，
# Y^2XY^2，Y^3。在poly（）输出的结果矩阵中，变量名称取自两个变量的
#阶数，例如，"1.2"表示X^1*Y^2。这里获得的是原始的多项式，如果构建正#交多项式，设定参数raw=FALSE
# 基于9个多项式项的RDA 
mite.trend.rda <- rda(mite.h ~ ., data=as.data.frame(mite.poly))
# 计算校正R2
(R2adj.poly <- RsquareAdj(mite.trend.rda)$adj.r.squared)
# 基于地理坐标正交的三阶项RDA
mite.poly.ortho <- poly(as.matrix(mite.xy), degree=3)
colnames(mite.poly.ortho) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")
mite.trend.rda.ortho <- rda(mite.h~., data=as.data.frame(mite.poly.ortho))
(R2adj.poly <- RsquareAdj(mite.trend.rda.ortho)$adj.r.squared)
# 使用Blanchet等（2008a）双终止准则的变量前向选择
(mite.trend.fwd <- forward.sel(mite.h, mite.poly.ortho, adjR2thresh=R2adj.poly))
# 新的RDA只保留6项 
(mite.trend.rda2 <- rda(mite.h ~ .,data=as.data.frame(mite.poly)[,mite.trend.fwd[,2]]))
# 全部典范轴检验和单轴检验
anova.cca(mite.trend.rda2, step=1000)
anova.cca(mite.trend.rda2, step=1000, by="axis")
# 绘制三个独立显著的空间模型（典范轴）图。以"sr.value"代替"s.value"生成图
# 中气泡图。
mite.trend.fit <- scores (mite.trend.rda2, choices=c(1,2,3), display="lc", scaling=1)
par(mfrow=c(1,3))
sr.value(mite.xy,mite.trend.fit[,1])
sr.value(mite.xy,mite.trend.fit[,2])
sr.value(mite.xy,mite.trend.fit[,3])
#如果将多项式项直接输入rda函数，代码可以这样写：mite.trend.rda <- #rda(mite.h~Xm+ym+I(Xm^2)+ I(Xm*Ym)+ I(Ym^2))；注意高阶项和混合项必
#须使用I（）函数区别，否则R会将这些项当作方差分析处理。
# 甲螨数据去趋势
# **************
anova(rda(mite.h, mite.xy))  # 结果：趋势显著
# 甲螨数据线性去趋势计算
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))





# 逐步构建PCNM变量
# ******************
# 1.一维取样：100个等距离取样点的样带，相邻两点距离为1
tr100 <- seq(1:100)            # 产生样线点
tr100.d1 <- dist(tr100)         #欧氏距离矩阵
thresh <- 1                    #设定削减阈值为1 
# 根据阈值削减欧氏距离矩阵
tr100.d1[tr100.d1 > thresh] <- 4*thresh  
# 削减距离矩阵的PCoA
tr100.PCoA <- cmdscale(tr100.d1, eig=TRUE, k=length(tr100)-1)
# 计算正特征根的数量 
(nb.ev <- length(which(tr100.PCoA$eig > 0.0000001)))
# PCNM变量矩阵
tr100.PCNM <- tr100.PCoA$points[,1:nb.ev]
# 绘制一些模拟正空间相关的PCNM变量(图7.3)
par(mfrow=c(4,2))
somePCNM <- c(1, 2, 4, 8, 15, 20, 30, 40)
for(i in 1:length(somePCNM)){
  plot(tr100.PCNM[,somePCNM[i]], type="l", ylab=c("PCNM", somePCNM[i]))
}
# 2.二维取样：等距离栅格。削减阈值设定为1，也可以设定为一个栅格4个点组
# 成的正方形对角线的距离，即sqrt(2) 
xygrid2 <- expand.grid(1:20, 1:20)
xygrid2.d1 <- dist(xygrid2)
thresh <- 1                    # 设定阈值距离为1
# 根据阈值削减欧氏距离矩阵
xygrid2.d1[xygrid2.d1>thresh] <- 4*thresh  
# 削减后距离矩阵的PCoA
xygrid2.PCoA <- cmdscale(xygrid2.d1, eig=TRUE, k=nrow(xygrid2)-1)
# 计算正特征根的数量   
(nb.ev2 <- length(which(xygrid2.PCoA$eig > 0.0000001)))
# PCNM变量矩阵
xygrid2.PCNM <- xygrid2.PCoA$points[,1:nb.ev2]
# 绘制一些模拟正空间相关的PCNM变量(图7.4)
par(mfrow=c(4,2))
somePCNM2 <- c(1, 2, 5, 10, 20, 50, 100, 150)
for(i in 1:length(somePCNM2)){
sr.value(xygrid2, xygrid2.PCNM[,somePCNM2[i]], method="greylevel", csize=0.35, sub=somePCNM2[i], csub=2)
}


# 甲螨数据的PCNM分析
# *******************
# 1a.逐步构建PCNM变量矩阵
# -------------------------------------
xy.d1 <- dist(mite.xy)
#搜索削减阈值：使用vegan程序包内spantree函数确定保证所有点都能够连接
#时的最小距离。使用等于或大于最小距离作为削减阈值。
spanning <- spantree(xy.d1)
dmin <- max(spanning$dist)
# 削减距离矩阵
xy.d1[xy.d1 > dmin] <- 4*dmin
# 削减后距离矩阵的PCoA
xy.PCoA <- cmdscale(xy.d1, k=nrow(mite.xy)-1, eig=TRUE)
# 计算正特征根的数量（PCNM可以模拟正或负的空间相关）
(nb.ev <- length(which(xy.PCoA$eig > 0.0000001)))
# 构建PCNM变量数据框
mite.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(mite.xy), 1:nb.ev])
# 1b. ...或自动构建PCNM变量
# ------------------------------------
# library(PCNM) # 如果还未加载PCNM程序包
xy.d1 <- dist(mite.xy)
mite.PCNM.auto <- PCNM(xy.d1)
summary(mite.PCNM.auto)
# 绘制确定削减阈值的最小拓展树
plot(mite.PCNM.auto$spanning, mite.xy)
(dmin <- mite.PCNM.auto$thresh) # 削减距离
(nb.ev <- length(mite.PCNM.auto$values)) # 特征根数量
# PCNM变量Moran指数（由第一距离等级0到削减阈值）；也见PCNM（）函数
# 产生的图（此处无显示图）
# Moran指数的期望值（代表无空间相关）
mite.PCNM.auto$expected_Moran
mite.PCNM.auto$Moran_I
# 正空间相关的特征函数
(select <- which(mite.PCNM.auto$Moran_I$Positive == TRUE))
length(select)  # I > E(I)条件下PCNM变量的数量
mite.PCNM.pos <- as.data.frame(mite.PCNM.auto$vectors)[,select]
# 1c. ... 或使用vegan程序包内pcnm（）函数
# --------------------------------------
mite.PCNM.vegan <- pcnm(dist(mite.xy))
mite.PCNM <- as.data.frame(mite.PCNM.vegan$vectors)
#dmin <- mite.PCNM.vegan$threshold
nb.ev <- length(which(mite.PCNM.vegan$values > 0.0000001))
#vegan包内pcnm（）函数获得的特征向量已经除以特征根的平方根，这样也
#可以作为空间变量。与PCNM（）函数不同，vegan包内pcnm（）函数不提供#Moran系数，如果需要选择仅仅保留正空间相关的特征函数，需要单独计算#Moran。
#系数。

# 2.运行基于去趋势甲螨数据的全模型PCNM分析
# ------------------------------------------------------------
mite.PCNM.rda <- rda(mite.h.det, mite.PCNM.pos)
anova.cca(mite.PCNM.rda)
# 3.如果分析为显著，计算校正R2和进行PCNM变量的前向选择
(mite.R2a <- RsquareAdj(mite.PCNM.rda)$adj.r.squared)
(mite.PCNM.fwd <- forward.sel(mite.h.det, as.matrix(mite.PCNM.pos), 
  adjR2thresh=mite.R2a))
#按照校正R2的原则，如果运行至PCNM5时获得的校正R2将略大于全模型的
#校正R2，但这个结果可能并不严谨。 
 (nb.sig.PCNM <- nrow(mite.PCNM.fwd)) # 显著PCNM的数量
# 依次排列显著的PCNM的变量
(PCNM.sign <- sort(mite.PCNM.fwd[,2]))
# 将显著的PCNM变量设定为新的对象
PCNM.red <- mite.PCNM.pos[,c(PCNM.sign)]
# 4.只用10个显著的PCNM变量进行新的PCNM分析（前向选择后R2adj=0.2713）
mite.PCNM.rda2 <- rda(mite.h.det ~ ., data=PCNM.red)
(mite.fwd.R2a <- RsquareAdj(mite.PCNM.rda2)$adj.r.squared)
anova.cca(mite.PCNM.rda2)
axes.test <- anova.cca(mite.PCNM.rda2, by="axis")
(nb.ax <- length(which(axes.test[,5] <= 0.05))) #显著轴的数量
# 5.绘制两个显著典范轴
mite.PCNM.axes <- scores (mite.PCNM.rda2, choices=c(1,2), display="lc", 
  scaling=1)
par(mfrow=c(1,2))
sr.value(mite.xy, mite.PCNM.axes[,1]) # ade4程序包函数：s.value
sr.value(mite.xy, mite.PCNM.axes[,2]) # ade4程序包函数：s.value
# 空间变化解读：两个显著典范轴与环境变量的回归分析
shapiro.test(resid(lm(mite.PCNM.axes[,1] ~ ., data=mite.env))) # 残差正态性检验
mite.PCNM.axis1.env <- lm(mite.PCNM.axes[,1]~., data=mite.env)
summary(mite.PCNM.axis1.env)
shapiro.test(resid(lm(mite.PCNM.axes[,2] ~ ., data=mite.env))) # 残差正态性检验
mite.PCNM.axis2.env <- lm(mite.PCNM.axes[,2] ~ ., data=mite.env)
summary(mite.PCNM.axis2.env)
#可以发现，与两个空间排序轴有关联的环境变量并不相同（除了灌丛shrubs
#这个变量）。当然，每个排序轴也不能被环境变量完全解释。更准确评估解
#释比例需要使用变差分解方法。
# 10个显著的PCNM变量图
# *************************
par(mfrow=c(2,5))
for(i in 1:ncol(PCNM.red)){
  sr.value(mite.xy, PCNM.red[,i], sub=PCNM.sign[i], csub=2)
}
# 甲螨PCNM分析-宽尺度
# *********************
mite.PCNM.broad <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(1,3,4)])
anova.cca(mite.PCNM.broad)
axes.broad <- anova.cca(mite.PCNM.broad, by="axis")
nb.ax.broad <- length(which(axes.broad[,5] <= 0.05))
nb.ax.broad                        # 显著轴的数量
# 绘制两个显著典范轴
mite.PCNMbroad.axes <- scores(mite.PCNM.broad, choices=c(1,2), display="lc", scaling=1)
par(mfrow=c(1,2))
s.value(mite.xy, mite.PCNMbroad.axes[,1])
s.value(mite.xy, mite.PCNMbroad.axes[,2])

# 解读空间变化：两个典范轴与环境变量的回归
mite.PCNMbroad.ax1.env <- lm(mite.PCNMbroad.axes[,1] ~ ., data=mite.env)
summary(mite.PCNMbroad.ax1.env)
mite.PCNMbroad.ax2.env <- lm(mite.PCNMbroad.axes[,2] ~ ., data=mite.env)
summary(mite.PCNMbroad.ax2.env)
#回归分析的结果清楚表明宽尺度空间结构与微地形和灌丛缺失密切相关。
# 甲螨PCNM分析-中尺度
# ********************
mite.PCNM.med <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(5,6,7,10,11)])
anova.cca(mite.PCNM.med)
axes.med <- anova.cca(mite.PCNM.med, by="axis")
(nb.ax.med <- length(which(axes.med[,5] <= 0.05)))  # 显著轴的数量
# 绘制显著典范轴（第二轴是边缘（marginally）显著）
mite.PCNMmed.axes <- scores(mite.PCNM.med, choices=c(1,2), display="lc", 
  scaling=1)
par(mfrow=c(1,2))
s.value(mite.xy, mite.PCNMmed.axes[,1])
s.value(mite.xy, mite.PCNMmed.axes[,2])
# 解读空间变化：两个典范轴与环境变量的回归
mite.PCNMmed.ax1.env <- lm(mite.PCNMmed.axes[,1] ~ ., data=mite.env)
summary(mite.PCNMmed.ax1.env)
mite.PCNMmed.ax2.env <- lm(mite.PCNMmed.axes[,2] ~ ., data=mite.env)
summary(mite.PCNMmed.ax2.env)
#回归分析的结果表明中尺度空间结构对应于土壤覆盖类型和土壤湿度（变量#WatrCont）
# 甲螨PCNM分析-微尺度
# ********************
mite.PCNM.fine <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(20,23)])
anova.cca(mite.PCNM.fine)
axes.fine <- anova.cca(mite.PCNM.fine, by="axis")
(nb.ax.fine <- length(which(axes.fine[,5] <= 0.05))) # 典范轴的数量
# 绘制显著典范轴
mite.PCNMfine.axes <- scores(mite.PCNM.fine, choices=1, display="lc", scaling=1)
par(mfrow=c(1,2))
s.value(mite.xy, mite.PCNMfine.axes)
# 解读空间变化：两个典范轴与环境变量的回归
mite.PCNMfine.ax1.env <- lm(mite.PCNMfine.axes ~ ., data=mite.env)
summary(mite.PCNMfine.ax1.env)
#回归分析的结果表明微尺度空间结构与环境变量有微弱的关系。仅有土壤
#覆盖类型（植物凋落物）是显著的。



# 使用quickPCNM（）函数进行一站式PCNM分析
# ****************************************
mite.PCNM.quick <- quickPCNM(mite.h, mite.xy)
summary(mite.PCNM.quick)
mite.PCNM.quick[[2]]   # 特征根
mite.PCNM.quick[[3]]   # 变量前向选择的结果
# quickPCNM（）函数并没有保留变量前向选择过程刚导致超过全模型的校
#正R2的PCNM5。因此，这里PCNM5并没有被包括在当前PCNM分析结果中。
#从quickPCNM函数输出结果，提取和绘制RDA结果（2型标尺）
# *****************************************************
plot(mite.PCNM.quick$RDA, scaling=2)
sp.scores1 <- scores(mite.PCNM.quick$RDA,  choices=1:2, scaling=2, display="sp")
arrows(0, 0, sp.scores1[,1], sp.scores1[,2], length=0, lty=1, col="red")
#2型标尺能够比较准确地展示某些物种与PCNM变量之间的关系。这些关系
#可以用于探索在哪个尺度上物种分布具有空间结构。




# 甲螨-趋势-环境-PCNM变差分解
# *****************************
# 1.检验趋势，如果显著，对坐标进行前向选择
# -----------------------------------------------------------
mite.XY.rda <- rda(mite.h, mite.xy)
anova.cca(mite.XY.rda)
(mite.XY.R2a <- RsquareAdj(mite.XY.rda)$adj.r.squared)
(mite.XY.fwd <- forward.sel(mite.h, as.matrix(mite.xy), 
  adjR2thresh=mite.XY.R2a))
XY.sign <- sort(mite.XY.fwd$order)
# 将显著的坐标变量赋予新的对象
XY.red <- mite.xy[,c(XY.sign)]
# 2. 环境变量检验和前向选择
# -------------------------------------
# 将环境变量3-5重新编码成二元变量
substrate <- model.matrix(~mite.env[,3])[,-1]
shrubs <- model.matrix(~mite.env[,4])[,-1]
topo <- model.matrix(~mite.env[,5])[,-1]
mite.env2 <- cbind(mite.env[,1:2], substrate, shrubs, topo)
# 环境变量的前向选择
mite.env.rda <- rda(mite.h, mite.env2)
(mite.env.R2a <- RsquareAdj(mite.env.rda)$adj.r.squared)
mite.env.fwd <- forward.sel(mite.h, mite.env2, adjR2thresh=mite.env.R2a,
 nperm=9999)
env.sign <- sort(mite.env.fwd$order)
env.red <- mite.env2[,c(env.sign)]
colnames(env.red)
# 3. PCNM变量的前向选择
# ----------------------------------
# 运行未去趋势甲螨数据的全模型PCNM分析
mite.undet.PCNM.rda <- rda(mite.h, mite.PCNM.pos)
anova.cca(mite.undet.PCNM.rda)
# 如果分析表明显著，计算校正R2和运行PCNM变量前向选择
(mite.undet.PCNM.R2a <- RsquareAdj(mite.undet.PCNM.rda)$adj.r.squared)
(mite.undet.PCNM.fwd <- forward.sel(mite.h, as.matrix(mite.PCNM.pos), 
  adjR2thresh=mite.undet.PCNM.R2a))
# 根据R2a准则，如果保留12个PCNM变量，获得的校正R2已经稍大于全模
# 型的校正R2。但这个"稍微超过"也是可行，并不一定很严格。
(nb.sig.PCNM <- nrow(mite.undet.PCNM.fwd)) # 显著的PCNM变量的数量
# 按顺序排列显著的PCNM变量
(PCNM.sign <- sort(mite.undet.PCNM.fwd$order))
# 赋予所有显著PCNM变量一个新的对象
PCNM.red <- mite.PCNM.pos[,c(PCNM.sign)]
# 4. 将显著PCNM变量任意分成宽尺度和微尺度变量
# ---------------------------------------------------------------------
# 宽尺度: PCNMs 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
PCNM.broad <- PCNM.red[,1:10]
# 微尺度: PCNMs 16, 20
PCNM.fine <- PCNM.red[,11:12]
# 5.甲螨-环境-趋势-PCNM变差分解
# -----------------------------------------------------------
(mite.varpart <- varpart(mite.h, env.red, XY.red, PCNM.broad, PCNM.fine))
par(mfrow=c(1,2))
showvarparts(4)
plot(mite.varpart, digits=2)
# 检验单独解释部分[a], [b], [c] 和 [d] 
# ********************************
# [a]部分，环境变量单独解释部分
anova.cca(rda(mite.h, env.red, cbind(XY.red, PCNM.broad, PCNM.fine)))
# [b]部分，趋势单独解释部分
anova.cca(rda(mite.h, XY.red, cbind(env.red, PCNM.broad, PCNM.fine))) 
# [c]部分，宽尺度空间变量单独解释部分
anova.cca(rda(mite.h, PCNM.broad, cbind(env.red, XY.red, PCNM.fine)))
# [d]部分，微尺度空间变量单独解释部分
anova.cca(rda(mite.h, PCNM.fine, cbind(env.red, XY.red, PCNM.broad)))
# 仅有环境变量和宽尺度空间变量单独解释部分显著


# 去趋势甲螨数据MEM分析
# **********************
# 最优空间权重矩阵的筛选
# **********************
# 1.基于Delaunay三角网搜索
# 使用mite.h.det作为响应变量，使用mite.del作为构建Delaunay三角网数据。
# 非权重矩阵（二元权重）；使用test.W函数从基于Delaunay三角网构建的MEM
# 变量中选择变量
# Delaunay 三角网 
(mite.del <- tri2nb(mite.xy)) 
mite.del.res <- test.W(mite.h.det, mite.del)
#屏幕结果显示最优模型的AICC值为-94.2（基于7个MEM变量）。
# 总结最优模型的输出结果
summary(mite.del.res$best)
# 最优模型未校正的R2
# 返回最小AICc值模型的R2
 (R2.del <- mite.del.res$best$R2[which.min(mite.del.res$best$AICc)])
# 最优模型的校正R2 (n = 70 和 m = 7)
RsquareAdj(R2.del,70,7)
# 2.通过距离函数设定权重的Delaunay三角网
# 距离归一化（除以最大值）并用幂参数alpha不同代表构建不同的初始W矩阵
f2 <- function(D, dmax, y) { 1 - (D/dmax)^y }
# 属于Delaunay三角网内连接的最大欧氏距离
max.d1 <- max(unlist(nbdists(mite.del, as.matrix(mite.xy)))) 
# 设定幂参数由2到10
mite.del.f2 <- test.W(mite.h.det, mite.del, f=f2, y=2:10, dmax=max.d1, 
  xy=as.matrix(mite.xy))
# 屏幕显示最优模型的AICC值为-95.4（基于6个MEM变量）
# 最优模型未校正的R2
(R2.delW <- mite.del.f2$best$R2[which.min(mite.del.f2$best$AICc)])
# 最优模型校正的R2 (n = 70 和 m = 6)
RsquareAdj(R2.delW,70,6)
# 3a.基于距离的连接矩阵（以点为中心、以距离阈值为半径范围内的点均连接）
# 通过去趋势甲螨数据的多元变异函数图（基于20个距离等级）选择距离阈值 
(mite.vario <- variogmultiv(mite.h.det, mite.xy, nclass=20))
plot(mite.vario$d, mite.vario$var, ty='b', pch=20, xlab="距离等级", ylab="C(距离)")
# 10个邻体矩阵（nb类型）构建
# 10个距离阈值向量
(thresh10 <- seq(give.thresh(dist(mite.xy)), 4, le=10))
# 生成10个邻体矩阵
# 每个邻体矩阵包含所有连接距离≤阈值
list10nb <- lapply(thresh10, dnearneigh, x=as.matrix(mite.xy), d1=0)
# 显示第一个邻体矩阵部分内容
print(listw2mat(nb2listw(list10nb[[1]], style="B"))[1:10,1:10], digits=1)
# 现在可以用test.W（）函数分析10个邻体矩阵，此处没有连接的权重
mite.thresh.res <- lapply(list10nb, test.W, Y=mite.h.det)
# 最小AICc值，最优模型，最优模型的距离阈值
mite.thresh.minAIC <- sapply(mite.thresh.res, function(x) min(x$best$AICc, 
  na.rm=TRUE))
min(mite.thresh.minAIC)  # 最小AICc (10个模型当中的最优模型)
which.min(mite.thresh.minAIC) # 最优模型的位置thresh10[which.min(mite.thresh.minAIC)]  # 削减距离阈值
#确认以最低AICc模型作为最优模型。这个模型的距离范围是多少？多少MEM 
#变量被选择？
# 3b.与上面一样属于第三类MEM模型，但连接通过距离幂函数1-(d/dmax)^y设
# 定连接权重矩阵
mite.thresh.f2 <- lapply(list10nb, function(x) test.W(x, Y=mite.h.det, f=f2, 
  y=2:10, dmax=max(unlist(nbdists(x, as.matrix(mite.xy)))), 
  xy=as.matrix(mite.xy))) 
# 最低AIC，最优模型
mite.f2.minAIC <- sapply(mite.thresh.f2, function(x) min(x$best$AICc, 
  na.rm=TRUE)) 
# 最小AICc (10个备选模型中的最优模型)
min(mite.f2.minAIC) 
# 最优模型的位置       
(nb.bestmod <- which.min(mite.f2.minAIC))  
# 最优模型实际的dmax 
(dmax.best <- mite.thresh.f2[nb.bestmod][[1]]$all[1,2]) 
# 最优MEM模型的信息提取
# ***********************
mite.MEM.champ <- unlist(mite.thresh.f2[which.min(mite.f2.minAIC)], 
  recursive=FALSE)
summary(mite.MEM.champ)
mite.MEM.champ$best$values   #特征根
mite.MEM.champ$best$ord      # 按照R2大小对MEM变量排位
# 最优模型保留的MEM变量
MEMid <- mite.MEM.champ$best$ord[1:which.min(mite.MEM.champ$best$AICc)]
sort(MEMid)
MEM.all <- mite.MEM.champ$best$vectors
MEM.select <- mite.MEM.champ$best$vectors[, sort(c(MEMid))]
colnames(MEM.select) <- sort(MEMid)
# 最优模型未校正的R2 
R2.MEMbest <- mite.MEM.champ$best$R2[which.min(mite.MEM.champ$best$AICc)]
# 最优模型校正的R2
RsquareAdj(R2.MEMbest, nrow(mite.h.det), length(MEMid))
# 使用plot.links（）函数绘制连接图
plot.links(mite.xy, thresh=dmax.best)
# 7个显著的MEM变量图 
# ********************
par(mfrow=c(2,4))
for(i in 1:ncol(MEM.select)){
  s.value(mite.xy,MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}
# 被7个MEM变量约束的去趋势甲螨数据RDA分析
# *****************************************
(mite.MEM.rda <- rda(mite.h.det~., as.data.frame(MEM.select)))
(mite.MEM.R2a <- RsquareAdj(mite.MEM.rda)$adj.r.squared)
anova.cca(mite.MEM.rda)
axes.MEM.test <- anova.cca(mite.MEM.rda, by="axis")
(nb.ax <- length(which(axes.MEM.test[,5] <= 0.05))) # Number of significant axes
# 绘制两个显著典范轴的样方坐标图
mite.MEM.axes <- scores (mite.MEM.rda, choices=c(1,2), display="lc", 
  scaling=1)
par(mfrow=c(1,2))
s.value(mite.xy, mite.MEM.axes[,1])
s.value(mite.xy, mite.MEM.axes[,2])
# 7个显著MEM变量的空间分布地图
# **************************
par(mfrow=c(2,4))
for(i in 1:ncol(MEM.select)){
  s.value(mite.xy, MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}
# MEM和PCNM变量相关性分析
# ----------------------------------------
cor(MEM.select, PCNM.red)
# 按照嵌套关系排列的连接矩阵
# Delaunay 三角网（与前面案例一样）
mite.del <- tri2nb(mite.xy) 
# Gabriel图
mite.gab <- graph2nb(gabrielneigh(as.matrix(mite.xy)), sym=TRUE)
# 相对邻体图（relative neighbourhood graph）
mite.rel <- graph2nb(relativeneigh(as.matrix(mite.xy)), sym=TRUE)
# 最小拓展树（minimum spanning tree）
mite.mst <- mst.nb(dist(mite.xy))	
# 邻体矩阵按照连接数排列
# 绘制连接矩阵图
par(mfrow=c(2,2))
plot(mite.del, mite.xy, col="red", pch=20, cex=1)
title(main="Delaunay三角网 ")
plot(mite.gab, mite.xy, col="purple", pch=20, cex=1)
title(main="Gabriel 图")
plot(mite.rel, mite.xy, col="dark green", pch=20, cex=1)
title(main="相对邻体")
plot(mite.mst, mite.xy, col="brown", pch=20, cex=1)
title(main="最小拓展树")
# 编辑连接
# ********
# 1.人机交互方式
plot(mite.del, mite.xy, col="red", pch=20, cex=2)
title(main="Delaunay triangulation")
mite.del2 <- edit.nb(mite.del, mite.xy)
# 如果要删除某个连接，点击连接的两个节点并按照屏幕提示操作。 
# 2.或者将nb类型对象转换为可编辑矩阵，然后通过命令代码剔除连接
mite.del.mat <- nb2mat(mite.del, style="B")
# 剔除样方23和样方35之间的连接
mite.del.mat[23,35] <- 0
mite.del.mat[35,23] <- 0
# 重新转换为nb对象
mite.del3 <- neig2nb(neig(mat01=mite.del.mat))
plot(mite.del3, mite.xy)
# 案例：Delaunay三角网内样方23的邻体列表
mite.del[[23]]   # 删除连接之前
mite.del2[[23]]  # 通过人机交互方法删除连接之后
mite.del3[[23]]  #通过命令代码方法删除连接之后
# 基于某个距离（给定半径）的连接矩阵
#使用PCNM案例中最小的削减距离dmin = 1.011187m作为半径 
mite.thresh4 <- dnearneigh(as.matrix(mite.xy), 0, dmin*4)
nb2mat(mite.thresh4)[1:10,1:10] # 显示矩阵的一部分
# 使用更短的距离 (1*dmin, 2*dmin)
mite.thresh1 <- dnearneigh(as.matrix(mite.xy), 0, dmin*1)
mite.thresh2 <- dnearneigh(as.matrix(mite.xy), 0, dmin*2)
#使用更长的距离
mite.thresh8 <- dnearneigh(as.matrix(mite.xy), 0, dmin*8)
#绘制部分连接矩阵
par(mfrow=c(1,2))
plot(mite.thresh1, mite.xy, col="red", pch=20, cex=0.8)
title(main="1 * dmin")
plot(mite.thresh4, mite.xy, col="red", pch=20, cex=0.8)
title(main="4 * dmin")
#1 * dmin版本的图显示有一个样方（样方7）没有被连接。为了避免出现这样
#的问题，可以使用略大一些的距离阈值。例如此处使用1.00111188可以使样
#方7也被连接。距离小于4m的样方对非常多，所以4 * dmin的图显得非常拥挤。
# "nb"对象转换为"listw"对象
# 案例：上面生成mite.thresh4对象。"B"代表二元矩阵
mite.thresh4.lw <- nb2listw(mite.thresh4, style="B")
print(listw2mat(mite.thresh4.lw)[1:10,1:10], digits=1)
# 生成空间权重矩阵W =矩阵B和矩阵A的Hadamard积
# ------------------------------------------------------------------
# 将连接矩阵内"1"以欧氏距离代替
mite.thresh4.d1 <- nbdists(mite.thresh4, as.matrix(mite.xy))
# 使用相对距离的反数（即1-相对距离）作为权重
mite.inv.dist <- lapply(mite.thresh4.d1, function(x) 1-x/max(dist(mite.xy)))
# 生成空间权重矩阵W。参数"B"代表二元数据但仅表达连接矩阵，而非权重
# 矩阵
mite.invdist.lw <- nb2listw(mite.thresh4, glist=mite.inv.dist, style="B")
print(listw2mat(mite.invdist.lw)[1:10,1:10], digits=2)
#使用相对欧氏距离的反数对所有样方对之间的非零连接进行权重赋值
# MEM变量计算（基于listw类型的对象）
# -----------------------------------------------
mite.invdist.MEM <- scores.listw(mite.invdist.lw, echo=TRUE)
summary(mite.invdist.MEM)
mite.invdist.MEM$values
barplot(mite.invdist.MEM$values)
# 每个特征向量Moran指数的检验
(mite.MEM.Moran <- test.scores(mite.invdist.MEM, mite.invdist.lw, 999))
# 显著空间相关的MEM变量
which(mite.MEM.Moran[,2] <= 0.05) 
length(which(mite.MEM.Moran[,2] <= 0.05))
#MEM变量31至35具有显著的Moran指数；但只有1，2，3，4，5，6，7（6
#和7的显著性依赖于置换的结果）具有正的空间相关
# 将MEM向量存储到新的对象
# 所有的MEM变量
mite.invdist.MEM.vec <- mite.invdist.MEM$vectors
# 空间正相关的MEM变量
MEM.Moran.pos <- which(mite.MEM.Moran[,1] > -1/(nrow(mite.invdist.MEM$vectors)-1))
mite.invdist.MEM.pos <- mite.invdist.MEM.vec[,MEM.Moran.pos]
# 显著的正空间相关MEM变量
MEM.Moran.pos.sig <- MEM.Moran.pos[which(mite.MEM.Moran[MEM.Moran.pos,2] <= 0.05)]
mite.invdist.MEM.pos.sig <- mite.invdist.MEM.vec[,MEM.Moran.pos.sig]
# MEM特征根vs. Moran指数散点图
plot(mite.invdist.MEM$values, mite.MEM.Moran$stat, ylab="Moran's I", 
  xlab="Eigenvalues")
text(-1, 0.5, paste("Correlation=", cor(mite.MEM.Moran$stat, 
  mite.invdist.MEM$values)))
# AEM分析
# *********
# 生成河流树形图代码 
# 参见Legendre和Legendre（1998，第47页）
lake1 <- c(1,0,1,1,0,0,0,0)
lake2 <- c(1,0,1,0,0,0,0,0)
lake3 <- c(1,1,0,0,0,0,0,0)
lake4 <- c(0,0,0,0,1,0,1,1)
lake5 <- c(0,0,0,0,0,1,1,1)
lake6 <- c(0,0,0,0,0,0,0,1)
arbor <- rbind(lake1, lake2, lake3, lake4, lake5, lake6)
# AEM变量构建
(arbor.aem <- aem(binary.mat=arbor))
arbor.aem.vec <- arbor.aem$vectors
# AEM特征函数也可以利用奇异值分解（函数svd（））获得，其实aem（）函数
#也做相同的运算。 
arbor.c = scale(arbor, center=TRUE, scale=FALSE)
arbor.svd = svd(arbor.c)
#奇异值
arbor.svd$d[1:5]    
# AEM特征函数
arbor.svd$u[,1:5]   
# 取样设计代码：10个跨河的样带，每个样带有4个样方，边缘的权重等比于相
# 对距离的平方的反数
# X-Y坐标
xy <- cbind(1:40, expand.grid(1:4, 1:10))
# nb类型对象（spdep程序包）包含类似国际象棋的"女王"的连接
nb <- cell2nb(4, 10, "queen")
# 样方-边缘矩阵（产生一个虚拟对象"0"）
edge.mat <- build.binary(nb, xy)
# 欧氏距离矩阵
D1.mat <- as.matrix(dist(xy))
# 提取边缘，剔除那些直接连接"0"的边缘
edges.b <- edge.mat$edges[-1:-4,]
# 构建一个代表每个边缘长度的向量
length.edge <- vector(length=nrow(edges.b))
for(i in 1:nrow(edges.b)){
	length.edge[i] <- D1.mat[edges.b[i,1], edges.b[i,2]]
}
# 将相对距离平方的反数作为每个边缘的权重
weight.vec <- 1-(length.edge/max(length.edge))^2
# 从edge.mat对象构建AEM特征函数 
example.AEM <- aem(build.binary=edge.mat, weight=weight.vec, rm.link0=TRUE)
example.AEM$values
ex.AEM.vec <- example.AEM$vectors
# 构建5种虚拟物种
# 两种随机分布的物种
sp12 <- matrix(trunc(rnorm(80,5,2),0),40)
# 其中一个物种限制分布于溪流上半部分
sp3 <- c(trunc(rnorm(20,8,2.5),0), rep(0,20))
# 另一个物种限制分布与样带的左半部
sp4 <- t(matrix(c(trunc(rnorm(20,8,3),0), rep(0,20)),10))
sp4 <- c(sp4[,1], sp4[,2], sp4[,3], sp4[,4], sp4[,5], sp4[,6], sp4[,7], 
  sp4[,8], sp4[,9], sp4[,10])
# 还有一个物种限制分布于左上部4个样方
sp5 <- c(4,7,0,0,3,8, rep(0,34))
sp <- cbind(sp12, sp3, sp4, sp5)
# 前20个AEM变量的全模型AEM分析（计算R2a）
#******************************************
AEM.20 <- rda(sp ~ ., as.data.frame(ex.AEM.vec[,1:20]))
R2a.AEM <- RsquareAdj(AEM.20)$adj.r.squared
AEM.fwd <- forward.sel(sp,ex.AEM.vec, adjR2thresh=R2a.AEM)
(AEM.sign <- sort(AEM.fwd[,2])) 
#将显著的AEM变量赋予新的对象
AEM.sign.vec <- ex.AEM.vec[,c(AEM.sign)] 
(sp.AEMsign.rda <- rda(sp ~ ., data=as.data.frame(AEM.sign.vec))) #显著AEM变
#量的RDA
anova.cca(sp.AEMsign.rda)
AEM.rda.axes.test <- anova.cca(sp.AEMsign.rda, by="axis")
#显著轴的数量
(nb.ax.AEM <- length(which(AEM.rda.axes.test[,5] <= 0.05))) 
# 绘制显著典范轴
AEM.rda.axes <- scores(sp.AEMsign.rda, choices=c(1,2), display="lc", 
  scaling=1)
par(mfrow=c(1,nb.ax.AEM))
for(i in 1:nb.ax.AEM) s.value(xy[,c(2,3)], AEM.rda.axes[,i])


# 多尺度排序（MSO）
# ****************
# 未去趋势的甲螨数据的MSO vs. 环境RDA
# ---------------------------------------------------
mite.undet.env.rda <- rda(mite.h, mite.env2)
mite.env.rda.mso <- mso(mite.undet.env.rda, mite.xy, grain=dmin, perm=999)
msoplot(mite.env.rda.mso, alpha=0.05/7)
mite.env.rda.mso
# 未去趋势甲螨数据MSO分析vs.RDA分析（MEM变量控制空间结构）
# ----------------------------------------------------------------------------------------
mite.undet.env.MEM <- rda(mite.h, mite.env2, as.data.frame(MEM.select))
mite.env.MEM.mso <- mso(mite.undet.env.MEM, mite.xy, grain=dmin, perm=999)
msoplot(mite.env.MEM.mso, alpha=0.05/7)
mite.env.MEM.mso
# 基于去趋势甲螨数据和环境数据的MSO分析
# ------------------------------------------------------------
# 对甲螨数据Y方向去趋势
mite.h.det2 <- resid(lm(as.matrix(mite.h) ~ mite.xy[,2]))
#对环境数据Y方向去趋势
env2.det <- resid(lm(as.matrix(mite.env2) ~ mite.xy[,2]))
# RDA 和 MSO
mitedet.envdet.rda <- rda(mite.h.det2, env2.det)
miteenvdet.rda.mso <- mso(mitedet.envdet.rda, mite.xy, grain=dmin, perm=999)
msoplot(miteenvdet.rda.mso, alpha=0.05/7)
miteenvdet.rda.mso
# 去趋势甲螨数据MSO分析vs.RDA分析（MEM控制空间结构）
# ----------------------------------------------------------------------
mite.det.env.MEM <- rda(mite.h.det2, env2.det, as.data.frame(MEM.select))
mite.env.MEM.mso <- mso(mite.det.env.MEM, mite.xy, grain=dmin, perm=999)
msoplot(mite.env.MEM.mso, alpha=0.05/7)
mite.env.MEM.mso
