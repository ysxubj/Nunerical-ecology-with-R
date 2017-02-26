# 导入CSV数据文件
# *****************
# 导入物种（群落）数据框（鱼类多度数据）
spe <- read.csv("DoubsSpe.csv", row.names=1)
# 导入环境数据框
env <- read.csv("DoubsEnv.csv", row.names=1)
# 导入空间坐标数据框
spa <- read.csv("DoubsSpa.csv", row.names=1)
# 基础函数
# ********
spe	#在控制台显示整个数据框的内容，但对于大样本的数据框
                 #并不建议直接显示
spe[1:5,1:10]	   #只展示前5行和前10列
head(spe)			#只展示前几行
nrow(spe)			#提取数据框总行数
ncol(spe)			#提取数据框总列数
dim(spe)		    #提取数据框的维度（显示数据框多少行，多少列）
colnames(spe)		#提取列名，在这里是物种名
rownames(spe)	    #提取行名，在这里一行代表一个样方
summary(spe)		 #以列为单位，对列变量进行描述性统计
#比较多度的中位值和平均值。大部分是对称分布吗？
# 多度数据总体分布情况
# *******************
# 整个多度数据值的范围
range(spe)
# 计算每种多度值的数量
ab <- table(unlist(spe))
ab
# 所有种混和在一起的多度分布柱状图
barplot(ab, las=1, xlab="多度等级", ylab="频度", col=gray(5:0/5))
# 多度数据中0值的数量
sum(spe==0)
# 多度数据中0值所占比例
sum(spe==0) / (nrow(spe)*ncol(spe))
#请观察多度频率分布柱状图，如何解读为什么0值（缺失）在数据框内频
#率这么高?
# 样方位置地图
# **************
# 生成空的绘图窗口（横纵坐标轴比例1:1，带标题）
# 从spa数据框获取地理坐标x和y
plot(spa, asp=1, type="n", main="样方位置",
	xlab="x坐标 (km)", ylab="y坐标 (km)")
# 加一条连接各个样方点的蓝色线（代表Doubs河）
lines(spa, col="light blue")
# 添加每个样方的编号
text(spa, row.names(spa), cex=0.8, col="red")
# 添加文本
text(70, 10, "上游", cex=1.2, col="red")
text(20, 120, "下游", cex=1.2, col="red")
#某些鱼类的分布地图
# ******************
# 将绘图窗口分割为4个绘图区域，每行两个
par(mfrow=c(2,2))
plot(spa, asp=1, col="brown", cex=spe$TRU, main="褐鳟",
	xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$OMB, main="茴鱼",
	xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$BAR, main="鲃鱼",
	xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$BCO, main="欧鳊",
xlab="x坐标 (km)", ylab="y坐标 (km)")
lines(spa, col="light blue")
#观察所生成的4张图，你就会明白为什么Verneaux 选择这4种鱼类作为不     #同区域的生态指示种，看了后面将要展示的环境因子空间分布情况会更清楚。
# 比较物种频度
# **************
# 计算每个物种出现的样方数
#按列进行计数，因此函数apply()第二个参数MARGIN应该设定为2
spe.pres <- apply(spe > 0, 2, sum)
# 按照升序的方式重新排列结果
sort(spe.pres)
# 计算频度百分比
spe.relf <- 100*spe.pres/nrow(spe)
round(sort(spe.relf), 1)	# 设置排列结果为1位小数
#绘柱状图
par(mfrow=c(1,2))		# 将绘图窗口垂直一分为二
hist(spe.pres, main="物种出现数", right=FALSE, las=1,
	xlab="出现数", ylab="物种数量",
	breaks=seq(0,30,by=5), col="bisque")
hist(spe.relf, main="物种相对频度", right=FALSE, las=1,
	xlab="出现率(%)", ylab="物种数量",
  	breaks=seq(0, 100, by=10), col="bisque")
 # 样方比较：物种丰富度
# ********************
# 计算每个样方内物种数
# 以行汇总，apply()函数第二个参数MARGIN应该设定为1
sit.pres <- apply(spe > 0, 1, sum)
#按照升序的方式重新排列结果
sort(sit.pres)
par(mfrow=c(1,2))	#将绘图窗口垂直一分为二
# 绘制样方沿着河流的分布位置和所含物种丰富度
plot(sit.pres,type="s", las=1, col="gray",
	main="物种丰富度-上下游的梯度",
	xlab="样方沿着河流的位置", ylab="物种丰富度")
text(sit.pres, row.names(spe), cex=.8, col="red")
# 使用地理坐标绘制气泡地图
plot(spa, asp=1, main="物种丰富度地图", pch=21, col="white",
	bg="brown", cex=5*sit.pres/max(sit.pres), xlab="x坐标 (km)",
	ylab="y坐标 (km)")
lines(spa, col="light blue")
#你能否辨析沿着河流哪里是物种丰富度的热点地区？
#计算生物多样性指数
# *****************
# 载入所需要的vegan程序包
library(vegan) # 如果未载入，需要执行这一步
#访问diversity() 帮助界面
?diversity
N0 <- rowSums(spe > 0)               #物种丰富度
H <- diversity(spe)                    # Shannon熵指数
N1 <- exp(H)                        # Shannon 多样性指数
N2 <- diversity(spe, "inv")              # Simpson多样性指数
J <- H/log(N0)                          # Pielou 均匀度
E1 <- N1/N0                            # Shannon均匀度 (Hill比率)
E2 <- N2/N0                            # Simpson均匀度 (Hill比率)
div <- data.frame(N0, H, N1, N2, E1, E2, J)
div
# 数据转化和标准化
##################
#访问decostand()帮助文件
?decostand
# 简单转化
# **********************
# 显示原始数据某一部分（多度数据）
spe[1:5, 2:4]
# 将多度数据转化为有-无（1-0）数据
spe.pa <- decostand(spe, method="pa")
spe.pa[1:5, 2:4]
#物种水平：两个方法；
#有-无数据或多度数据
# *******************
# 通过每个数值除以该物种最大值标准化多度
# 注意: 这里参数MARGIN=2 (默认值)
spe.scal <- decostand(spe, "max")
spe.scal[1:5,2:4]
# 计算每列最大值
apply(spe.scal, 2, max)
#这些标准化过程是否正确运行？最好利用绘图函数或总结函数密切追踪.。
#通过每个数值除以该物种总和标准化多度（每个物种的相对多度）
#注意: 这里需要设定参数MARGIN=2
spe.relsp <- decostand(spe, "total", MARGIN=2)
spe.relsp[1:5,2:4]
#计算标准化后数据每列总和
apply(spe.relsp, 2, sum)
# 样方水平：3种方法；有-无数据或多度数据
# ***************************************
#通过每个数值除以该样方总和标准化多度 (每个样方相对多度或相对频度)
#注意: 这里参数MARGIN=1 (默认值)
spe.rel <- decostand(spe, "total")	# 默认MARGIN = 1
spe.rel[1:5,2:4]
#计算标准化后数据每列总和以检验标准化的过程是否正确
apply(spe.rel, 1, sum)
#赋予每个行向量长度（范数）为1（即平方和为1）.
spe.norm <- decostand(spe, "normalize")
spe.norm[1:5,2:4]
# 验证每个行向量的范数
norm <- function(x) sqrt(x%*%x)
apply(spe.norm, 1, norm)
#这个转化也称为"弦转化"：如果用欧氏距离函数去计算弦转化后的数据，#将获得弦距离矩阵（见第3章）。在PCA和RDA（见第5、6章）及k-means
#聚类（见第4章）分析前通常需要对数据进行弦转化。
# 计算相对频度（样方层面），然后取平方根
spe.hel <- decostand(spe, "hellinger")
spe.hel[1:5,2:4]
# 计算标准化后数据每行向量的范数
apply(spe.hel,1,norm)
#这个转化也称为Hellinger转化。如果用欧氏距离函数去计算Hellinger转
#化后的数据，将获得Hellinger距离矩阵（见第3章）。在PCA和RDA（见
#第5、6章）及k-means聚类（见第4章）分析前通常需要对数据进行Hellinger
#转化。注意，Hellinger转化等同于数据先平方根转化后再进行弦转化。
# 物种和样方同时标准化
# ****************************
# 卡方转化
spe.chi <- decostand(spe, "chi.square")
spe.chi[1:5,2:4]
# 请查看没有物种的样方8转化后将会怎样？
spe.chi[7:9,]
#如果用欧氏距离函数去计算卡方转化后的数据，将获得卡方距离矩阵（见
#第3章）
# Wisconsin标准化：多度数据首先除以该物种最大值后再除以该样方总和
spe.wis <- wisconsin(spe)
spe.wis[1:5,2:4]
# 常见种（石泥鳅 stone loach）转化后的多度箱线图
# *******************************************
par(mfrow=c(2,2))
boxplot(spe$LOC, sqrt(spe$LOC), log1p(spe$LOC),
	las=1, main="简单转化",
	names=c("原始数据", "sqrt", "log"), col="bisque")
boxplot(spe.scal$LOC, spe.relsp$LOC,
	las=1, main="物种标准化",
	names=c("max", "total"), col="lightgreen")
boxplot(spe.hel$LOC, spe.rel$LOC, spe.norm$LOC,
	las=1, main="样方标准化",
	names=c("Hellinger", "total", "norm"), col="lightblue")
boxplot(spe.chi$LOC, spe.wis$LOC,
	las=1, main="双标准化",
	names=c("Chi-square", "Wisconsin"), col="orange")
#比较多度数据转化或标准化前后的数据分布范围和分布情况。
#绘制物种从河流上游到下游分布图
# ******************************
par(mfrow=c(2,2))
plot(env$das, spe$TRU, type="l", col=4, main="Raw data",
  xlab="Distance from the source [km]", ylab="Raw abundance code")
lines(env$das, spe$OMB, col=3)
lines(env$das, spe$BAR, col="orange")
lines(env$das, spe$BCO, col=2)
lines(env$das, spe$LOC, col=1, lty="dotted")
plot(env$das, spe.scal$TRU, type="l", col=4, main="Species profiles (max)",
  xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.scal$OMB, col=3)
lines(env$das, spe.scal$BAR, col="orange")
lines(env$das, spe.scal$BCO, col=2)
lines(env$das, spe.scal$LOC, col=1, lty="dotted")

plot(env$das, spe.hel$TRU, type="l", col=4,
  main="Site profiles (Hellinger)",
  xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.hel$OMB, col=3)
lines(env$das, spe.hel$BAR, col="orange")
lines(env$das, spe.hel$BCO, col=2)
lines(env$das, spe.hel$LOC, col=1, lty="dotted")

plot(env$das, spe.chi$TRU, type="l", col=4,
  main="Double profiles (Chi-square)",
  xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.chi$OMB, col=3)
lines(env$das, spe.chi$BAR, col="orange")
lines(env$das, spe.chi$BCO, col=2)
lines(env$das, spe.chi$LOC, col=1, lty="dotted")
legend("topright", c("Brown trout", "Grayling", "Barbel", "Common bream",
  "Stone loach"), col=c(4,3,"orange",2,1), lty=c(rep(1,4),3))
#比较这些图，并解释它们的不同。
# 部分环境变量的气泡地图
# *******************************************
par(mfrow=c(2,2))
plot(spa, asp=1, main="海拔", pch=21, col="white", bg="red",
	cex=5*env$alt/max(env$alt), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="流量", pch=21, col="white", bg="blue",
	cex=5*env$deb/max(env$deb), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="氧含量", pch=21, col="white", bg="green3",
	cex=5*env$oxy/max(env$oxy), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="硝酸盐浓度", pch=21, col="white", bg="brown",
	cex=5*env$nit/max(env$nit), xlab="x", ylab="y")
lines(spa, col="light blue")
 #哪幅图最能展示上下游的梯度？如何解释其他环境变量的空间分布格局？
#线条图
# *****
par(mfrow=c(2,2))
plot(env$das, env$alt, type="l", xlab="离源头距离 (km)",
	ylab="海拔 (m)", col="red", main="海拔")
plot(env$das, env$deb, type="l", xlab="离源头距离 (km)",
	ylab="流量 (m3/s)", col="blue", main="流量")
plot(env$das, env$oxy, type="l", xlab="离源头距离 (km)",
	ylab="氧含量 (mg/L)", col="green3", main="氧含量")
plot(env$das, env$nit, type="l", xlab="离源头距离 (km)",
	ylab="硝酸盐浓度 (mg/L)", col="brown", main="硝酸盐浓度")
# 所有变量对之间的二维散点图
# **************************
#载入自编的函数R脚本
source("panelutils.R")  # panelutils.R脚本文件必须与当前R工作空间在同一文件
#夹下
# 带频度分布的柱状图和光滑拟合曲线的双变量散点图
op <- par(mfrow=c(1,1), pty="s")
pairs(env, panel=panel.smooth, diag.panel=panel.hist,
	main="双变量散点图（带频度分布图和平滑曲线）")
par(op)
 #从柱状图能否看出哪些变量符合正态分布？
#需要注意的是，对于回归分析和典范排序，并没有要求解释变量符合正态
#分布。是否有很多散点图显示出变量之间的线性关系或至少是单调关系？
# 某个环境变量简单转化
# ********************
range(env$pen)
# 坡度变量对数转化(y = ln(x))
# 比较转化前后数值的柱状图和箱线图
par(mfrow=c(2,2))
hist(env$pen, col="bisque", right=FALSE, main="坡度频度分布图",ylab="频度",xlab="坡度")
hist(log(env$pen), col="light green", right=F, main="对数化坡度频度分布图",ylab="频度",xlab="对数化坡度")
boxplot(env$pen, col="bisque", main="坡度箱线图", ylab="坡度")
boxplot(log(env$pen), col="light green", main="对数化坡度箱线图",
	ylab="对数化坡度")
# 所有环境变量的标准化
# *********************
# 中心化和标准化=标准化变量 (z-scores)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)	# 平均值 = 0
apply(env.z, 2, sd)		# 标准差 = 1
# 使用scale()函数也可以运行相同的标准化(输出的是矩阵)
env.z <- as.data.frame(scale(env))
