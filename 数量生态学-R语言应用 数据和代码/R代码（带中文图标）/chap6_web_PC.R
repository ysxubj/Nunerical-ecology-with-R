# CHAPTER 6 - CANONICAL ORDINATION
# ********************************
# 载入所需程序包
library(ade4)
library(vegan)
library(packfor)
# 此程序包可以从 https://r-forge.r-project.org/R/?group_id=195 下载
# 如果是MacOS X系统，packfor程序包内forward.sel函数的运行需要加载
# gfortran程序包。用户必须从"cran.r-project.org"网站内选择"MacOS X"，
# 然后选择"tools"安装gfortran程序包。
library(MASS)
library(ellipse)
library(FactoMineR)
# 附加函数
source("evplot.R")
source("hcoplot.R")
# 导入CSV数据文件
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# 删除没有数据的样方8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]
# 提取环境变量das（离源头距离）以备用
das <- env[, 1]
# 从环境变量矩阵剔除das变量
env <- env[, -1]
# 将slope变量（pen）转化为因子（定性）变量
pen2 <- rep("very_steep", nrow(env))
pen2[env$pen <= quantile(env$pen)[4]] = "steep"
pen2[env$pen <= quantile(env$pen)[3]] = "moderate"
pen2[env$pen <= quantile(env$pen)[2]] = "low"
pen2 <- factor(pen2, levels=c("low", "moderate", "steep", "very_steep"))
table(pen2)
# 生成一个含定性坡度变量的环境变量数据框env2
env2 <- env
env2$pen <- pen2
# 将所有解释变量分为两个解释变量子集
# 地形变量（上下游梯度）子集
envtopo <- env[, c(1:3)]
names(envtopo)
#水体化学属性变量子集
envchem <- env[, c(4:10)]
names(envchem)
# 物种数据Hellinger转化
spe.hel <- decostand(spe, "hellinger")
# 基于Hellinger 转化的鱼类数据RDA，解释变量为对象env2包括的环境变量
# 关注省略模式的公式
spe.rda <- rda(spe.hel~., env2) 
summary(spe.rda) # 2型标尺（默认）
#这里使用一些默认的选项，即 scale=FALSE（基于协方差矩阵的RDA）和#scaling=2
#如何从rda（）输出结果中获得典范系数
coef(spe.rda)
# 提取校正R2
# **********
# 从rda的结果中提取未校正R2 
(R2 <- RsquareAdj(spe.rda)$r.squared)
# 从rda的结果中提取校正R2
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)
#可以看出，校正R2总是小于R2。校正R2作为被解释方差比例的无偏估计，后#面的变差分解部分所用的也是校正R2。
# RDA三序图
# **********
# 1型标尺：距离三序图
plot(spe.rda, scaling=1, main="RDA三序图：spe.hel～env2 - 1型标尺- 加权和样方坐标")
#此排序图同时显示所有的元素：样方、物种、定量解释变量（用箭头表示）
#和因子变量的形心。为了与定量解释变量区分，物种用不带箭头的线表示。
spe.sc <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")
# 2型标尺（默认）：相关三序图
plot(spe.rda, main="RDA三序图：spe.hel～env2 - 2型标尺- 加权和样方坐标")
spe2.sc <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe2.sc[, 1], spe2.sc[, 2], length=0, lty=1, col="red")
# 样方坐标是环境因子线性组合 
# 1型标尺
plot(spe.rda, scaling=1, display=c("sp", "lc", "cn"), 
	main="RDA三序图：spe.hel～env2 - 1型标尺- 拟合的样方坐标")
arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=0, lty=1, col="red")
# 2型标尺
plot(spe.rda, display=c("sp", "lc", "cn"), 
	main="RDA三序图：spe.hel～env2 - 2型标尺- 拟合的样方坐标")
arrows(0, 0, spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col="red")
# RDA所有轴置换检验
anova.cca(spe.rda, step=1000)	
# 每个典范轴逐一检验
anova.cca(spe.rda, by="axis", step=1000)
#每个典范轴的检验只能输入由公式模式获得的rda结果。有多少个轴结果是
#显著的？
# 使用Kaiser-Guttman准则确定残差轴
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]
#很明显，还有一部分有意思的变差尚未被目前所用的这套环境变量解释。


# 偏RDA：固定地形变量影响后，水体化学属性的效应
# 简单模式：X和W可以是分离的定量变量表格
spechem.physio <- rda(spe.hel, envchem, envtopo)
spechem.physio
summary(spechem.physio)
# 公式模式；X和W必须在同一数据框内
class(env)
spechem.physio2 <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo 
 + Condition(alt + pen + deb), data=env)
spechem.physio2
#上面两个分析的结果完全相同。
# 偏RDA检验（使用公式模式获得的RDA结果，以便检验每个轴）
anova.cca(spechem.physio2, step=1000)
anova.cca(spechem.physio2, step=1000, by="axis")
# 偏RDA三序图（使用拟合值的样方坐标）
# 1型标尺
plot(spechem.physio, scaling=1, display=c("sp", "lc", "cn"), 
main="RDA三序图：spe.hel～chem ︳Tope- 1型标尺- 拟合的样方坐标")
spe3.sc <- scores(spechem.physio, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe3.sc[, 1], spe3.sc[, 2], length=0, lty=1, col="red")
# 2型标尺
plot(spechem.physio, display=c("sp", "lc", "cn"), 
main="RDA三序图：spe.hel～chem ︳Tope- 2型标尺- 拟合的样方坐标")
spe4.sc <- scores(spechem.physio, choices=1:2, display="sp")
arrows(0, 0, spe4.sc[,1], spe4.sc[,2], length=0, lty=1, col="red")



# 两个RDA结果的变量方差膨胀因子(VIF)
# ********************************************
# 本章第一个RDA结果：包括所有环境因子变量
vif.cca(spe.rda)
vif.cca(spechem.physio) # 偏RDA
# 使用双终止准则（Blanchet等，2008a）前向选择解释变量
# 1.包含所有解释变量的RDA全模型 
spe.rda.all <- rda(spe.hel ~ ., data=env)
# 2.全模型校正R2
(R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)
# 3.使用packfors 包内forward.sel（）函数选择变量
# library(packfor) #如果尚未载入packfors包，需要运行这一步
forward.sel(spe.hel, env, adjR2thresh=R2a.all)
#注意，正如这个函数的提示信息所示，选择最后一个变量其实违背了
#adjR2thresh终止准则，所以pen变量最终不应该在被选变量列表中。
# 使用vegan包内ordistep（）函数前向选择解释变量。该函数可以对因子变量进# 行选择，也可以运行解释变量的逐步选择和后向选择。
step.forward <- ordistep(rda(spe.hel ~ 1, data=env), scope = formula(spe.rda.all ), 
direction="forward", pstep = 1000)
RsquareAdj(rda(spe.hel ~ alt, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt+oxy, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt+oxy+dbo, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt+oxy+dbo+pen, data=env))$adj.r.squared
#简约的RDA分析
# **************
spe.rda.pars <- rda(spe.hel ~ alt + oxy + dbo, data=env)
spe.rda.pars
anova.cca(spe.rda.pars, step=1000)
anova.cca(spe.rda.pars, step=1000, by="axis")
vif.cca(spe.rda.pars)
(R2a.pars <- RsquareAdj(spe.rda.pars)$adj.r.squared)
# 简约RDA的三序图（拟合的样方坐标）
# **********************************
# 1型标尺
plot(spe.rda.pars, scaling=1, display=c("sp", "lc", "cn"), 
	main="RDA三序图：spe.hel～alt+oxy+dbo - 1型标尺 - 拟合的样方坐标")
spe4.sc = scores(spe.rda.pars, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe4.sc[, 1], spe4.sc[, 2], length=0, lty=1, col="red")
# 2型标尺
plot(spe.rda.pars, display=c("sp", "lc", "cn"), 
	main="RDA三序图：spe.hel～alt+oxy+dbo - 2型标尺 - 拟合的样方坐标")
spe5.sc = scores(spe.rda.pars, choices=1:2, display="sp")
arrows(0, 0, spe5.sc[,1], spe5.sc[,2], length=0, lty=1, col="red")
#如果第三典范轴也显著，可以选择绘制轴1和轴3、轴2和轴3的三序图。
# 两组解释变量的变差分解
# ***********************
# 变差分解说明图
par(mfrow=c(1,3))
showvarparts(2) # 两组解释变量
showvarparts(3) #三组解释变量
showvarparts(4) #四组解释变量
# 1.带所有环境变量的变差分解 
spe.part.all <- varpart(spe.hel, envchem, envtopo)
spe.part.all
plot(spe.part.all, digits=2) 
#这些图内校正R2是正确的数字，但是韦恩图圆圈大小相同，未与R2的大小
#成比例。
# 2.分别对两组环境变量进行前向选择
spe.chem <- rda(spe.hel, envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
forward.sel(spe.hel, envchem, adjR2thresh=R2a.all.chem, nperm=9999)
spe.topo <- rda(spe.hel, envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
forward.sel(spe.hel, envtopo, adjR2thresh=R2a.all.topo, nperm=9999)
# 解释变量简约组合（基于变量选择的结果）
names(env)
envchem.pars <- envchem[, c(4,6,7)]
envtopo.pars <- envtopo[, c(1,2)]
# 变差分解
(spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars))
plot(spe.part, digits=2)
# 所有可测部分的置换检验
# [a+b]部分的检验
anova.cca(rda(spe.hel, envchem.pars), step=1000)
# [b+c]部分的检验
anova.cca(rda(spe.hel, envtopo.pars), step=1000)
# [a+b+c]部分的检验
env.pars <- cbind(envchem.pars, envtopo.pars)
anova.cca(rda(spe.hel, env.pars), step=1000)
# [a]部分的检验
anova.cca(rda(spe.hel, envchem.pars, envtopo.pars), step=1000)
# [c]部分的检验
anova.cca(rda(spe.hel, envtopo.pars, envchem.pars), step=1000)
#各个部分置换检验有不显著的吗？
# 3.无变量nit（硝酸盐浓度）的变差分解
envchem.pars2 <- envchem[, c(6,7)]
(spe.part2 <- varpart(spe.hel, envchem.pars2, envtopo.pars))
plot(spe.part2, digits=2)


# 基于RDA的双因素MANOVA
# **************************
# 生成代表海拔的因子变量（3个水平，每个水平含9个样方）
alt.fac <- gl(3, 9)
# 生成近似模拟pH值的因子变量
pH.fac <- as.factor(c(1, 2, 3, 2, 3, 1, 3, 2, 1, 2, 1, 3, 3, 2, 1, 1, 2, 3, 
  2, 1, 2, 3, 2, 1, 1, 3, 3))
# 两个因子是否平衡？
table(alt.fac, pH.fac)
# 用Helmert对照法编码因子和它们的交互作用项
alt.pH.helm <- model.matrix(~ alt.fac * pH.fac, 
 contrasts=list(alt.fac="contr.helmert", pH.fac="contr.helmert"))
alt.pH.helm
#检查当前对照法产生的表格，哪一列代表海拔因子、pH值和交互作用项？
# 检查Helmert对照表属性1：每个变量的和为1
apply(alt.pH.helm[, 2:9], 2, sum) 
# 检查Helmert对照表属性2：变量之间不相关
cor(alt.pH.helm[, 2:9]) 
# 使用函数betadisper（）（vegan包）（Marti Anderson检验）验证组内协方差矩阵# 的齐性
spe.hel.d1 <- dist(spe.hel[1:27,])
# 海拔因子
(spe.hel.alt.MHV <- betadisper(spe.hel.d1, alt.fac))
anova(spe.hel.alt.MHV)
permutest(spe.hel.alt.MHV) # 置换检验
# pH值因子
(spe.hel.pH.MHV <- betadisper(spe.hel.d1, pH.fac))
anova(spe.hel.pH.MHV)
permutest(spe.hel.pH.MHV) #置换检验
#组内协方差齐性，可以继续分析。
# 首先检验交互作用项。海拔因子和pH因子构成协变量矩阵 
interaction.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 6:9], alt.pH.helm[, 2:5]) 
anova(interaction.rda, step=1000, perm.max=1000) 
#交互作用是否显著？显著的交互作用表示一个因子的影响依赖另一个因子
#的水平，这将妨害主因子变量的分析。
# 检验海拔因子的效应，此时pH值因子和交互作用项作为协变量矩阵 
factor.alt.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 2:3], alt.pH.helm[, 4:9]) 
anova(factor.alt.rda, step=1000, perm.max=1000, strata=pH.fac) 
#海拔因子影响是否显著？
#检验pH值因子的效应，此时海拔值因子和交互作用项作为协变量矩阵
factor.pH.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 4:5], 
  alt.pH.helm[, c(2:3, 6:9)]) 
anova(factor.pH.rda, step=1000, perm.max=1000, strata=alt.fac) 
#pH值影响是否显著？
# RDA和显著影响的海拔因子三序图
alt.rda.out <- rda(spe.hel[1:27,]~., as.data.frame(alt.fac))
plot(alt.rda.out, scaling=1, display=c("sp", "wa", "cn"), 
	main="Multivariate ANOVA, factor altitude - scaling 1 - wa scores")
spe.manova.sc <- scores(alt.rda.out, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.manova.sc[, 1], spe.manova.sc[, 2], length=0, col="red")




# 基于距离的RDA分析（db-RDA）
# ****************************
# 1.分步计算
spe.bray <- vegdist(spe[1:27, ], "bray")
spe.pcoa <- cmdscale(spe.bray, k=nrow(spe[1:27, ])-1, eig=TRUE, add=TRUE)
spe.scores <- spe.pcoa$points
# 交互作用的检验。从协变量矩阵获得Helmert对照编码海拔因子和pH值因子
interact.dbrda <- rda(spe.scores[1:27, ], alt.pH.helm[, 6:9], alt.pH.helm[, 2:5]) 
anova(interact.dbrda, step=1000, perm.max=1000) 
#交互作用是否显著？如果没有，可以继续检验主因子的效应（此处未显示）
# 2.直接使用vegan包内capscale（）函数运行。只能以模型界面运行。响应变量
#可以是原始数据矩阵。
interact.dbrda2 <- capscale(spe[1:27,] ~ alt.fac*pH.fac + Condition(alt.fac+pH.fac), distance="bray", add=TRUE)
anova(interact.dbrda2, step=1000, perm.max=1000)
# 或者响应变量可以是相异矩阵。
interact.dbrda3 <- capscale(spe.bray ~ alt.fac*pH.fac + Condition(alt.fac+pH.fac), add=TRUE)
anova(interact.dbrda3, step=1000, perm.max=1000)

# 二阶解释变量的RDA
# *******************
# 生成das和das正交二阶项（由poly（）函数获得）矩阵
das.df <- poly(das, 2)
colnames(das.df) <- c("das", "das2")
# 验证两个变量是否显著
forward.sel(spe, das.df)
# RDA和置换检验
spe.das.rda <- rda(spe ~ ., as.data.frame(das.df))
anova(spe.das.rda)
# 三序图（拟合的样方坐标，2型标尺）
plot(spe.das.rda, scaling=2, display=c("sp", "lc", "cn"), 
	main="RDA三序图：spe～das + das2 - 2型标尺 - 拟合的样方坐标")
spe6.sc <- scores(spe.das.rda, choices=1:2, scaling=2, display="sp")
arrows(0, 0, spe6.sc[, 1], spe6.sc[, 2], length=0, lty=1, col="red")
# 4种鱼类的分布地图 
# ******************
par(mfrow=c(2, 2))
plot(spa$x, spa$y, asp=1, col="brown", cex=spe$TRU, 
	xlab="x (km)", ylab="y (km)", main="褐鳟")
lines(spa$x, spa$y, col="light blue")
plot(spa$x, spa$y, asp=1, col="brown", cex=spe$OMB, 
	xlab="x (km)", ylab="y (km)", main="鳟鱼")
lines(spa$x, spa$y, col="light blue")
plot(spa$x, spa$y, asp=1, col="brown", cex=spe$ABL, 
	xlab="x (km)", ylab="y (km)", main="欧鮊鱼")
lines(spa$x, spa$y, col="light blue")
plot(spa$x, spa$y, asp=1, col="brown", cex=spe$TAN, 
	xlab="x (km)", ylab="y (km)", main="鲤鱼")
lines(spa$x, spa$y, col="light blue")
自写代码角#4
为了能够正确自写RDA分析代码，有必要参考Legendre和Legendre 
（1998）第11.1节相关内容。下面是计算步骤（基于协方差矩阵的RDA）
   1.计算中心化的物种数据矩阵与标准化解释变量矩阵的多元线性回归，获得拟合值矩阵；
   2.计算拟合值矩阵的PCA；
   3.计算两类样方坐标；
   4.结果输出。
   下面代码解释部分使用的公式编码与Legendre和Legendre（1998）一致。
   下面的代码集中在RDA约束部分，目的是使读者对RDA数学过程感兴趣，而不是最优化程序。
myRDA <- function(Y,X) {
	# 1.数据的准备
	# *************
	Y.mat <- as.matrix(Y)
	Yc <- scale(Y.mat, scale=FALSE)
	X.mat <- as.matrix(X)
	Xcr  <- scale(X.mat)
	# 2.多元线性回归的计算
	# *********************
	# 回归系数矩阵 (eq. 11.4)
	B <- solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% Yc
	# 拟合值矩阵(eq. 11.5)
	Yhat <- Xcr %*% B
	# 残差矩阵
	Yres <- Yc - Yhat
	#维度
	n <- nrow(Y)
	p <- ncol(Y)
	m <- ncol(X)
	# 3. 拟合值PCA分析 
	# ******************
	# 协方差矩阵 (eq. 11.7)
	S <- cov(Yhat)
	# 特征根分解
	eigenS <- eigen(S)
	# 多少个典范轴？
	kc <- length(which(eigenS$values > 0.00000001))
	# 典范轴特征根
	ev <- eigenS$values[1:kc]
	# 矩阵Yc（中心化）的总方差（惯量）
	trace = sum(diag(cov(Yc)))
	# 正交特征向量（响应变量的贡献，1型标尺）
	U <- eigenS$vectors[,1:kc]
	row.names(U) <- colnames(Y)
	# 样方坐标（vegan包内'wa' 坐标，1型标尺eq. 11.12)
	F <- Yc %*% U
	row.names(F) <- row.names(Y)
	# 样方约束（vegan包内'lc' 坐标，2型标尺eq. 11.13)
	Z <- Yhat %*% U
	row.names(Z) <- row.names(Y)
	# 典范系数 (eq. 11.14)
	CC <- B %*% U
	row.names(CC) <- colnames(X)
	# 解释变量
	# 物种-环境相关
	corXZ <- cor(X,Z)
	# 权重矩阵的诊断
	D <- diag(sqrt(ev/trace))
	# 解释变量双序图坐标
	coordX <- corXZ %*% D    #1型标尺 
	coordX2 <- corXZ         #2型标尺
	row.names(coordX) <- colnames(X)
	row.names(coordX2) <- colnames(X)
	# 相对特征根平方根转化（为2型标尺）
	U2 <- U %*% diag(sqrt(ev))
	row.names(U2) <- colnames(Y)
	F2 <- F %*% diag(1/sqrt(ev))
	row.names(F2) <- row.names(Y)
	Z2 <- Z %*% diag(1/sqrt(ev))
	row.names(Z2) <- row.names(Y)
	# 未校正R2 
	R2 <- sum(ev/trace)
	# 校正R2
	R2adj <- 1-((n-1)/(n-m-1))*(1-R2)
	# 4.残差的PCA 
	# *******************
	# 与第5章相同，写自己的代码，可以从这里开始... 
	#     eigenSres <- eigen(cov(Yres))
	#     evr <- eigenSres$values
	# 5.输出Output
	result <- list(trace, R2, R2a,ev,CC,U,F,Z,coordX, U2,F2,Z2,coordX2)
	names(result) <- c("Total_variance", "R2", "R2a","Can_ev", "Can_coeff", 
    "Species_sc1", "wa_sc1", "lc_sc1", "Biplot_sc1", "Species_sc2", 
    "wa_sc2", "lc_sc2", "Biplot_sc2") 

	result
}
#将此函数应用到Doubs鱼类数据和环境数据的RDA分析
doubs.myRDA <- myRDA(spe.hel,env)
summary(doubs.myRDA)



# 典范对应分析（CCA）
# *******************
# 原始鱼类数据的CCA分析，解释变量为env2中所有环境变量
spe.cca <- cca(spe ~ ., env2)
spe.cca
summary(spe.cca) #2型标尺（默认）
# CCA三序图（使用拟合的样方坐标）
# ***********************************
# 1型标尺：物种坐标等比例于相对特征根 
# 样方坐标是物种坐标的加权平均
plot(spe.cca, scaling=1, display=c("sp","lc","cn"), main="CCA三序图：spe～env2 - 1型标尺")
# 2型标尺（默认）：样方坐标等比例于相对特征根  
# 物种坐标是样方坐标的加权平均
plot(spe.cca, display=c("sp","lc","cn"), main="CCA三序图：spe～env2 - 2型标尺")
# CCA无物种的双序图（1型标尺）(拟合的样方坐标)
# ***********************************************************
plot(spe.cca, scaling=1, display=c("lc", "cn"), 
 main="CCA双序图：spe～env2 - 1型标尺")
# CCA无样方的双序图（2型标尺）
# **********************************
plot(spe.cca, scaling=2, display=c("sp", "cn"), 
 main="CCA双序图：spe～env2 - 2型标尺")
# CCA结果的置换检验
# *******************
# 全模型置换检验
anova(spe.cca, step=1000)
# 每轴置换检验
anova(spe.cca, by="axis", step=1000)
# 使用vegan包ordistep（）函数对CCA模型的解释变量进行筛选#****************************************************
# ordistep（）函数允许使用因子变量，例如env2中的"pen"因子变量
cca.step.forward <- ordistep(cca(spe ~ 1, data=env2), scope=formula(spe.cca), 
direction="forward", pstep=1000)
# 仅使用alt、oxy和dbo三个解释变量的简约CCA分析
# *******************************************
(spe.cca.pars <- cca(spe ~ alt + oxy + dbo, data=env2))
anova.cca(spe.cca.pars, step=1000)
anova.cca(spe.cca.pars, step=1000, by="axis")
vif.cca(spe.cca)
vif.cca(spe.cca.pars)
# 可人机交互的三维排序图
# ***********************
# 样方排序图（wa坐标）
# ***********************
ordirgl(spe.cca.pars, type="t", scaling=1)
#使用鼠标的右键和滑轮可以缩放三维图的大小，使用左键可以转动三维图。
# 将加权平均的坐标连接到线性组合坐标
orglspider(spe.cca.pars, scaling=1, col="purple")
#紫色的连接线显示CCA模型的拟合数据的好坏。连接线越短，拟合越好。
# 附带聚类分析结果的样方（wa坐标）三维排序图
# *******************************************
# 样方编号不同颜色代表不同的聚类簇
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward"), 4)
ordirgl(spe.cca.pars, type="t", scaling=1, ax.col="black", col=gr+1)
# 将样方连接到聚类簇的形心
orglspider(spe.cca.pars, gr, scaling=1)
#样方沿着主要生态梯度很好被聚类。需要牢记当前是基于鱼类数据的分析。
# 完整的三维CCA三序图
# ***********************
ordirgl(spe.cca.pars, type="t", scaling=2)
orgltext(spe.cca.pars, display="species", type="t", scaling=2, col="cyan")
# 绘制物种组（使用Jaccard相似系数进行分组）
# ***********************************************************
gs <- cutree(hclust(vegdist(t(spe), method="jaccard"), "ward"), 4)
ordirgl(spe.cca.pars, display="species", type="t", col=gs+1)
rgl.quit（）  #关闭rgl工作系统


 # 线性判别式分析（LDA）
# *******************
# 基于Hellinger转化鱼类数据的样方Ward聚类分析（分4组）
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward"), 4)  #如果之前没有保存gr
env.pars2 <- as.matrix(env[, c(1, 9, 10)])
# 使用函数betadisper（）（vegan包）（Marti Anderson检验）验证解释变量组内# 协方差矩阵的齐性
env.pars2.d1 <- dist(env.pars2)
(env.MHV <- betadisper(env.pars2.d1, gr))
anova(env.MHV)
permutest(env.MHV) #置换检验
#组内协方差矩阵不齐，需要对alt和dbo两个变量进行对数化
env.pars3 <- cbind(log(env$alt), env$oxy, log(env$dbo))
colnames(env.pars3) <- c("alt.ln", "oxy", "dbo.ln") 
row.names(env.pars3) <- row.names(env)
env.pars3.d1 <- dist(env.pars3)
(env.MHV2 <- betadisper(env.pars3.d1, gr))
permutest(env.MHV2)
#组内协方差矩阵显示齐性，可以继续分析
# LDA计算（判别函数）
env.pars3.df <- as.data.frame(env.pars3)
(spe.lda <- lda(gr ~ alt.ln + oxy + dbo.ln, data=env.pars3.df))
# 此结果对象内包含大量解读LDA的必要信息
summary(spe.lda)
# 显示三个变量分组平均值
spe.lda$means
#计算范数标准化特征向量（matrix C，eq.11.33），即标准化判别函数系数
 (Cs <- spe.lda$scaling)
# 计算典范特征根
spe.lda$svd^2
# 计算在典范变量空间内样方的坐标
(Fp <- predict(spe.lda)$x)
# 替代方式: Fp <- scale(env.pars3.df, center=TRUE, scale=FALSE) %*% C
# 对象的分类
(spe.class <- predict(spe.lda)$class)
# 对象属于分类组的后验概率
 (spe.post <- predict(spe.lda)$posterior)
# 预设分类和预测分类的比较表格
(spe.table <- table(gr, spe.class))
# 正确分类的比例
diag(prop.table(spe.table, 1))
# 绘制典范变量空间内样方位置图，样方颜色代表不同的组
plot(Fp[, 1], Fp[, 2], type="n")
text(Fp[, 1], Fp[, 2], row.names(env), col=c(as.numeric(spe.class)+1))
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
# 绘制95%的椭圆图
for(i in 1:length(levels(as.factor(gr)))){
  cov <- cov(Fp[gr==i, ])
  centre <- apply(Fp[gr==i, ], 2, mean)
lines(ellipse(cov, centre= centre, level=0.95))
}
# 新样方的分类（识别） 
# 生成一个新的样方（ln(alt)=6.8, oxygen=90 和 ln(dbo)=3.2）
newo <- c(6.8, 90, 3.2)
newo <- as.data.frame(t(newo)) # 必须在同一行
colnames(newo) <- colnames(env.pars3)
(predict.new <- predict(spe.lda, newdata=newo))
# 基于刀切法（jackknife）分类的LDA（即留一法交叉验证)
spe.lda.jac <- lda(gr ~ alt.ln + oxy + dbo.ln, data=env.pars3.df, CV=TRUE)
summary(spe.lda.jac)
# 正确分类的数量和比例
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
diag(prop.table(spe.jac.table, 1))

# 典范相关分析（CCorA）
# **********************
# 数据的准备（对数据进行转化使变量分布近似对称）
envchem2 <- envchem
envchem2$pho <- log(envchem$pho)
envchem2$nit <- sqrt(envchem$nit)
envchem2$amm <- log1p(envchem$amm)
envchem2$dbo <- log(envchem$dbo)
envtopo2 <- envtopo
envtopo2$alt <- log(envtopo$alt)
envtopo2$pen <- log(envtopo$pen)
envtopo2$deb <- sqrt(envtopo$deb)
# CCorA (基于标准化的变量)
chem.topo.CCorA <- CCorA(envchem2, envtopo2, stand.Y=TRUE, stand.X=TRUE, 
  nperm=999)
chem.topo.CCorA
biplot(chem.topo.CCorA, plot.type="biplot")

# 协惯量分析（CoIA）
# *******************
# 两个环境变量矩阵的PCA排序 
dudi.chem <- dudi.pca(envchem2, scale=TRUE, scan=FALSE, nf=3)
dudi.topo <- dudi.pca(envtopo2, scale=TRUE, scan=FALSE, nf=2)
dudi.chem$eig/sum(dudi.chem$eig) # 每轴特征根比例dudi.topo$eig/sum(dudi.topo$eig) # 每轴特征根比例
all.equal(dudi.chem$lw,dudi.topo$lw) #两个分析每行权重是否相等？ 
# 协惯量分析
coia.chem.topo <- coinertia(dudi.chem,dudi.topo, scan=FALSE, nf=2)
coia.chem.topo
coia.chem.topo$eig[1]/sum(coia.chem.topo$eig) # 第1个特征根解释量
summary(coia.chem.topo)
randtest(coia.chem.topo, nrepet=999)    # 置换检验
plot(coia.chem.topo)

# 多元因子分析（MFA）
# *******************
# 三组变量的MFA
# 组合三个表格（Hellinger转化的物种数据、地形变量和水体化学属性）
tab3 <- data.frame(spe.hel, envtopo, envchem)
dim(tab3)
(grn <- c(ncol(spe), ncol(envtopo), ncol(envchem)))
# 计算MFA（附带多图）
# graphics.off（）  # 关闭前面的绘图窗口
t3.mfa <- MFA(tab3, group=grn, type=c("c","s","s"), ncp=2,
	name.group=c("鱼类群落","地形","水质"))
t3.mfa
plot(t3.mfa, choix="ind", habillage="none")
plot(t3.mfa, choix="ind", habillage="none", partial="all")
plot(t3.mfa, choix="var", habillage="group")
plot(t3.mfa, choix="axes")
# RV系数及检验
 (rvp <- t3.mfa$group$RV)
rvp[1,2] <- coeffRV(spe.hel, scale(envtopo))$p.value
rvp[1,3] <- coeffRV(spe.hel, scale(envchem))$p.value
rvp[2,3] <- coeffRV(scale(envtopo), scale(envchem))$p.value
round(rvp[-4,-4], 6)
# 特征根和方差百分百
t3.mfa$eig
ev <- t3.mfa$eig[,1]
names(ev) <- 1:nrow(t3.mfa$eig)
windows(title="MFA eigenvalues")
evplot(ev)
# 选择最典型的变量
aa <- dimdesc(t3.mfa, axes=1:2, proba=0.0001)
# 保留最显著（相关性）的变量排序图 
varsig <- t3.mfa$quanti.var$cor[unique(c(rownames(aa$Dim.1$quanti),
  rownames(aa$Dim.2$quanti))),]
plot(varsig[,1:2], asp=1, type="n", xlim=c(-1,1), ylim=c(-1,1))
abline(h=0, lty=3)
abline(v=0, lty=3)
symbols(0, 0, circles=1, inches=FALSE, add=TRUE)
arrows(0, 0, varsig[,1], varsig[,2], length=0.08, angle=20)
for (v in 1:nrow(varsig)) {
  if (abs(varsig[v,1]) > abs(varsig[v,2])) {
    if (varsig[v,1] >= 0) pos <- 4
    else pos <- 2
    }
  else {
    if (varsig[v,2] >= 0) pos <- 3
    else pos <- 1
    }
  text(varsig[v,1], varsig[v,2], labels=rownames(varsig)[v], pos=pos)
}        
