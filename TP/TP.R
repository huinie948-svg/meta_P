library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
library(orchaRd)
library(patchwork)
library(ggplot2)
library(dplyr) #加载dplyr包
library(ggpmisc) #加载ggpmisc包
library(ggpubr)
ap <- read.csv("TP.csv",fileEncoding = "latin1")
# 查看数据结构，确认各列的数据类型
str(ap)


# 将必要的列转换为数值型
# 使用as.numeric()转换，并处理可能的转换问题
ap$xc <- as.numeric(as.character(ap$xc))
ap$sc <- as.numeric(as.character(ap$sc))
ap$nc <- as.numeric(as.character(ap$nc))
ap$xe <- as.numeric(as.character(ap$xe))
ap$se <- as.numeric(as.character(ap$se))
ap$ne <- as.numeric(as.character(ap$ne))

# 检查转换后的结果
str(ap)
ap_effects <-  escalc(
  measure = "ROM",  # Ratio of Means (对数响应比)
  m1i = xe,         # 处理组均值
  m2i = xc,         # 对照组均值
  sd1i = se,  # 处理组标准差（已转换）
  sd2i = sc,  # 对照组标准差（已转换）
  n1i = ne,         # 处理组样本量
  n2i = nc,         # 对照组样本量
  data = ap
)
# 重命名效应量列
ap_effects <- ap_effects %>%
  rename(
    yi_lnRR = yi,      # natural-log response ratio
    vi_lnRR = vi       # sampling variance of lnRR
  )
# Backward-compatible aliases used by the legacy plotting code.
# In this project RR stores lnRR; it is not exponentiated.
ap_effects$RR <- ap_effects$yi_lnRR
ap_effects$Vi <- ap_effects$vi_lnRR
# 查看结果
head(ap_effects)
write.csv(ap_effects, "ap_effects.csv")


ap <- read.csv("ap_effects.csv", fileEncoding = "latin1")
# Check data
head(ap)

# 1. The number of Obversation
total_number <- nrow(ap)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 187 

# 2. The number of Study
unique_studyid_number <- length(unique(ap$study.id))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# umber of unique StudyID: 49 


#### 3. Overall effect size
total_effect_model <- rma.mv(
  yi = RR,
  V  = Vi,
  random = ~ 1 | study.id,
  data = ap,
  method = "REML",
  control = list(
    optimizer = "optim",
    optmethod = "BFGS",     # 也可以试 "Nelder-Mead"
    maxit     = 10000,
    rel.tol   = 1e-8
  )
)

summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016880.8629   4033761.7257   4033765.7257   4033772.1772   4033765.7913   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0343  0.1851     49     no  study.id 
# 
# Test for Heterogeneity:
#   Q(df = 186) = 4619884.2425, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub    
# 0.0182  0.0280  0.6498  0.5158  -0.0367  0.0732    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = ap, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# 
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# 
# Test for Funnel Plot Asymmetry: z = 0.3036, p = 0.7615
# Limit Estimate (as sei -> 0):   b = 0.0031 (CI: -0.0449, 0.0512)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# 
# Average Effect Size:         0.0076
# Amount of Heterogeneity:     0.0631
# Observed Significance Level: 0.7008
# Target Significance Level:   0.05
# 
# Fail-safe N: 0


p1 <- orchard_plot(total_effect_model,
                   mod = "1",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p1





#### 8. Subgroup analysis
### 8.1 methods
ap_method <- subset(ap, methods %in% c("acid digestion", "microwave digestion"))
#
ap_method$methods <- droplevels(factor(ap_method$methods))
# The number of Observations and StudyID
group_summary <- ap_method %>%
  group_by(methods) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
#A tibble: 2 × 3
# methods             Observations Unique_StudyID
# <fct>                      <int>          <int>
#   1 acid digestion                47             14
# 2 microwave digestion            8              1
overall_model_ap_method <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + methods, random = ~ 1 | study.id, data = ap_method, method = "REML")
# QM and p value
summary(overall_model_ap_method)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -249.9453   499.8906   505.8906   511.8015   506.3804   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0098  0.0992     15     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 53) = 927.0286, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1.9362, p-val = 0.3798
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# methodsacid digestion         0.0166  0.0304   0.5453  0.5856  -0.0430  0.0761    
# methodsmicrowave digestion   -0.1883  0.1471  -1.2802  0.2005  -0.4766  0.1000    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p2 <- orchard_plot(overall_model_ap_method,
                   mod = "methods",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p2
coef_rotation <- coef(overall_model_ap_method)
vcov_rotation <- vcov(overall_model_ap_method)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# methodsacid digestion methodsmicrowave digestion 
# "a"                        "a" 













### 8.2 Inoculation.location
ap_Inoculation.location <- subset(ap, inoculation.location %in% c("soils","seeds","seedings","roots"))
#
ap_Inoculation.location$inoculation.location <- droplevels(factor(ap_Inoculation.location$inoculation.location))
# The number of Observations and StudyID
group_summary <- ap_Inoculation.location %>%
  group_by(inoculation.location) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 4 × 3
# inoculation.location Observations Unique_StudyID
# <fct>                       <int>          <int>
#   1 roots                          10              2
# 2 seedings                        5              2
# 3 seeds                          41              8
# 4 soils                         121             38

overall_model_ap_Inoculation.location <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculation.location, random = ~ 1 | study.id, data = ap_Inoculation.location, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculation.location)
# Multivariate Meta-Analysis Model (k = 177; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016893.1934   4033786.3867   4033796.3867   4033812.1532   4033796.7460   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0292  0.1710     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 173) = 4529364.1978, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 5.6235, p-val = 0.2291
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# inoculation.locationroots       0.0950  0.1211   0.7844  0.4328  -0.1424  0.3324    
# inoculation.locationseedings   -0.1555  0.1333  -1.1664  0.2435  -0.4167  0.1058    
# inoculation.locationseeds      -0.0104  0.0326  -0.3186  0.7500  -0.0743  0.0535    
# inoculation.locationsoils       0.0282  0.0276   1.0210  0.3072  -0.0259  0.0823    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p3 <- orchard_plot(overall_model_ap_Inoculation.location,
                   mod = "inoculation.location",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p3
coef_rotation <- coef(overall_model_ap_Inoculation.location)
vcov_rotation <- vcov(overall_model_ap_Inoculation.location)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# inoculation.locationroots    inoculation.locationsoils    inoculation.locationseeds inoculation.locationseedings 
# "a"                          "a"                          "a"                          "a" 









### 8.3 Inoculant.type
ap_inoculant.type <- subset(ap, inoculant.type %in% c("bacteria", "fungi", "amf","mix"))
#
ap_inoculant.type$inoculant.type <- droplevels(factor(ap_Inoculant.type$inoculant.type))
# The number of Observations and StudyID
group_summary <- ap_inoculant.type %>%
  group_by(inoculant.type) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# A tibble: 4 × 3
# inoculant.type Observations Unique_StudyID
# <fct>                 <int>          <int>
#   1 amf                     329             78
# 2 bacteria                578            119
# 3 fungi                   116             30
# 4 mix                      91             20

overall_model_ap_Inoculant.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculant.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculant.type)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016880.4190   4033760.8381   4033770.8381   4033786.8855   4033771.1770   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0372  0.1928     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 183) = 4544301.0411, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 7.7734, p-val = 0.1002
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# inoculant.typeamf         0.0037  0.0362   0.1008  0.9197  -0.0674  0.0747    
# inoculant.typebacteria   -0.0256  0.0337  -0.7601  0.4472  -0.0916  0.0404    
# inoculant.typefungi       0.1431  0.0623   2.2960  0.0217   0.0209  0.2653  * 
#   inoculant.typemix        -0.0322  0.0392  -0.8197  0.4124  -0.1091  0.0447    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Extract model coefficients and covariance matrix
p4 <- orchard_plot(overall_model_ap_Inoculant.type,
                   mod = "inoculant.type",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p4
coef_rotation <- coef(overall_model_ap_Inoculant.type)
vcov_rotation <- vcov(overall_model_ap_Inoculant.type)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# inoculant.typefungi      inoculant.typeamf inoculant.typebacteria      inoculant.typemix 
# "a"                   "a"                    "a"                    "a" 










### 8.4 Experimental.type
ap_Experimental.type <- subset(ap, Experimental.type %in% c("Greenhouse", "Field"))
#
ap_Experimental.type$Experimental.type <- droplevels(factor(ap_Experimental.type$Experimental.type))
# The number of Observations and StudyID
group_summary <- ap_Experimental.type %>%
  group_by(Experimental.type) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)
# A tibble: 2 × 3
# Experimental.type Observations Unique_StudyID
# <fct>                    <int>          <int>
# 1 Field                      227             51
# 2 Greenhouse                 885            156

overall_model_ap_Experimental.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + experimental.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Experimental.type)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016880.6158   4033761.2315   4033767.2315   4033776.8926   4033767.3641   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0333  0.1824     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 185) = 4584342.9882, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1.3339, p-val = 0.5133
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# experimental.typefield        -0.0271  0.0552  -0.4910  0.6234  -0.1353  0.0811    
# experimental.typegreenhouse    0.0334  0.0320   1.0454  0.2959  -0.0292  0.0960    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Extract model coefficients and covariance matrix
p5 <- orchard_plot(overall_model_ap_Experimental.type,
                   mod = "experimental.type",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p5
coef_rotation <- coef(overall_model_ap_Experimental.type)
vcov_rotation <- vcov(overall_model_ap_Experimental.type)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# experimental.typegreenhouse      experimental.typefield 
# "a"                         "a" 











### 8.5 Inoculum
ap_Inoculum <- subset(ap, inoculum %in% c("bacteria yes","fungi yes", "amf yes", "amf inoculum yes", "mix yes"))
#
ap_Inoculum$inoculum <- droplevels(factor(ap_Inoculum$inoculum))
# The number of Observations and StudyID
group_summary <- ap_Inoculum %>%
  group_by(inoculum) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
#A tibble: 5 × 3
# inoculum         Observations Unique_StudyID
# <fct>                   <int>          <int>
#   1 amf inoculum yes            8              4
# 2 amf yes                    44             11
# 3 bacteria yes               96             28
# 4 fungi yes                  26              7
# 5 mix yes                     7              3

overall_model_ap_Inoculum <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculum, random = ~ 1 | study.id, data = ap_Inoculum, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculum)
# Multivariate Meta-Analysis Model (k = 181; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016886.3051   4033772.6101   4033784.6101   4033803.6330   4033785.1072   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0325  0.1803     47     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 176) = 4544443.6715, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 8.2379, p-val = 0.1436
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# inoculumamf inoculum yes    0.0046  0.0406   0.1142  0.9091  -0.0749  0.0842    
# inoculumamf yes             0.0097  0.0613   0.1587  0.8739  -0.1103  0.1298    
# inoculumbacteria yes       -0.0236  0.0350  -0.6754  0.4994  -0.0921  0.0449    
# inoculumfungi yes           0.1664  0.0690   2.4127  0.0158   0.0312  0.3015  * 
#   inoculummix yes            -0.0305  0.0410  -0.7430  0.4575  -0.1109  0.0499    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Extract model coefficients and covariance matrix
p6 <- orchard_plot(overall_model_ap_Inoculum,
                   mod = "inoculum",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p6
coef_rotation <- coef(overall_model_ap_Inoculum)
vcov_rotation <- vcov(overall_model_ap_Inoculum)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# inoculumfungi yes          inoculumamf yes inoculumamf inoculum yes     inoculumbacteria yes          inoculummix yes 
# "a"                     "a"                      "a"                      "a"                      "a" 










### 8.6 Inoculant.quantity
ap_Inoculant.quantity <- subset(ap, inoculant.quantity %in% c("yes", "no"))
#
ap_Inoculant.quantity$inoculant.quantity <- droplevels(factor(ap_Inoculant.quantity$inoculant.quantity))
# The number of Observations and StudyID
group_summary <- ap_Inoculant.quantity %>%
  group_by(inoculant.quantity) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 2 × 3
# inoculant.quantity Observations Unique_StudyID
# <fct>                     <int>          <int>
#   1 no                         1072            149
# 2 yes                         634             62

overall_model_ap_Inoculant.quantity <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculant.quantity, random = ~ 1 | study.id, data = ap_Inoculant.quantity, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculant.quantity)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -154803.1715   309606.3430   309612.3430   309628.6652   309612.3571   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0852  0.2919    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 620288.6887, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 14531.7729, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# inoculant.quantityno     0.3158  0.0220  14.3775  <.0001  0.2727  0.3588  *** 
#   inoculant.quantityyes    0.1302  0.0220   5.9203  <.0001  0.0871  0.1733  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ap_Inoculant.quantity)
vcov_rotation <- vcov(overall_model_ap_Inoculant.quantity)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# inoculant.quantityno inoculant.quantityyes 
# "a"                   "a" 










### 8.7 Soil.Origin
ap_Soil.Origin <- subset(ap, Soil.Origin %in% c("cropland", "forest", "Grassland", "plantation", "other"))
#
ap_Soil.Origin$Soil.Origin <- droplevels(factor(ap_Soil.Origin$Soil.Origin))
# The number of Observations and StudyID
group_summary <- ap_Soil.Origin %>%
  group_by(Soil.Origin) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)
# Soil.Origin Observations Unique_StudyID
# <fct>              <int>          <int>
#   1 cropland             277             58
# 2 forest                48              8
# 3 Grassland             18              1
# 4 other                689            117
# 5 plantation            34              8

overall_model_ap_Soil.origin <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + soil.origin, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Soil.origin)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016880.7689   4033761.5379   4033771.5379   4033787.5853   4033771.8768   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0364  0.1908     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 183) = 4578170.1157, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 1.0124, p-val = 0.9079
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# soil.origincropland     -0.0121  0.0727  -0.1668  0.8675  -0.1545  0.1303    
# soil.originforest        0.0803  0.0939   0.8551  0.3925  -0.1038  0.2643    
# soil.originother         0.0169  0.0349   0.4832  0.6289  -0.0516  0.0853    
# soil.originplantation    0.0160  0.1134   0.1412  0.8877  -0.2062  0.2383    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract model coefficients and covariance matrix
p7 <- orchard_plot(overall_model_ap_Soil.origin,
                   mod = "soil.origin",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p7
coef_rotation <- coef(overall_model_ap_Soil.origin)
vcov_rotation <- vcov(overall_model_ap_Soil.origin)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# soil.originplantation   soil.origincropland      soil.originother     soil.originforest  soil.origingrassland 
# "a"                   "a"                   "a"                   "a"                   "a"













### 8.9 extraction.part
ap_extraction.part <- subset(ap, extraction.part %in% c("rhizosphere", "bulk"))
#
ap_extraction.part$extraction.part <- droplevels(factor(ap_extraction.part$extraction.part))
# The number of Observations and StudyID
group_summary <- ap_extraction.part %>%
  group_by(extraction.part) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 2 × 3
# extraction.part Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 bulk                     108             27
# 2 rhizosphere               74             22

overall_model_ap_extraction.part <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + extraction.part, random = ~ 1 | study.id, data = ap_extraction.part, method = "REML")
# QM and p value
summary(overall_model_ap_extraction.part)
# Multivariate Meta-Analysis Model (k = 182; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016883.2528   4033766.5056   4033772.5056   4033782.0845   4033772.6420   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0353  0.1880     47     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 180) = 4091128.4966, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1.2909, p-val = 0.5244
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub    
# extraction.partbulk           0.0078  0.0317  0.2444  0.8069  -0.0544  0.0700    
# extraction.partrhizosphere    0.0348  0.0334  1.0430  0.2969  -0.0306  0.1002    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p9 <- orchard_plot(overall_model_ap_extraction.part,
                   mod = "extraction.part",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p9
coef_rotation <- coef(overall_model_ap_extraction.part)
vcov_rotation <- vcov(overall_model_ap_extraction.part)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# extraction.partrhizosphere        extraction.partbulk 
# "a"                        "a"








### 8.10 Sterilization
ap_Sterilization <- subset(ap, sterilization %in% c("sterilization", "unsterilization"))
#
ap_Sterilization$sterilization <- droplevels(factor(ap_Sterilization$sterilization))
# The number of Observations and StudyID
group_summary <-ap_Sterilization %>%
  group_by(sterilization) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# A tibble: 2 × 3
# sterilization   Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 sterilization             61             11
# 2 unsterilization          124             37

overall_model_ap_Sterilization <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + sterilization, random = ~ 1 | study.id, data = ap_Sterilization, method = "REML")
# QM and p value
summary(overall_model_ap_Sterilization)

# Multivariate Meta-Analysis Model (k = 185; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016882.5838   4033765.1676   4033771.1676   4033780.7960   4033771.3016   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0376  0.1940     48     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 183) = 4159448.3238, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 0.7906, p-val = 0.6735
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# sterilizationsterilization     -0.0109  0.0633  -0.1727  0.8629  -0.1351  0.1132    
# sterilizationunsterilization    0.0292  0.0335   0.8722  0.3831  -0.0365  0.0950    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p10 <- orchard_plot(overall_model_ap_Sterilization,
                   mod = "sterilization",          # 截距模型用 "1"
                   group = "study.id",
                   xlab = "Effect size(lnRR)",
                   trunk.size = 0.5,        # 组均值的置信区间粗细
                   branch.size = 2,         # 预测区间粗细
                   angle = 0,
                   legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p10
coef_rotation <- coef(overall_model_ap_Sterilization)
vcov_rotation <- vcov(overall_model_ap_Sterilization)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# sterilizationunsterilization   sterilizationsterilization 
# "a"                          "a"











### 8.10 Stress
ap_Stress <- subset(ap, Stress %in% c("Stress", "No Stress"))
#
ap_Stress$Stress <- droplevels(factor(ap_Stress$Stress))
# The number of Observations and StudyID
group_summary <-ap_Stress %>%
  group_by(Stress) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)
# # A tibble: 2 × 3
# Stress    Observations Unique_StudyID
# <fct>            <int>          <int>
#   1 No Stress          938            201
# 2 Stress             167             29

overall_model_ap_Stress <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + stress, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Stress)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -38855.7370   77711.4739   77717.4739   77727.1350   77717.6065   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0657  0.2564     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 185) = 2537930.3752, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 3956069.2639, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# stressno stress   -0.0310  0.0382  -0.8108  0.4175  -0.1058  0.0439      
# stressstress       0.5956  0.0382  15.5950  <.0001   0.5208  0.6705  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p11 <- orchard_plot(overall_model_ap_Stress,
                    mod = "stress",          # 截距模型用 "1"
                    group = "study.id",
                    xlab = "Effect size(lnRR)",
                    trunk.size = 0.5,        # 组均值的置信区间粗细
                    branch.size = 2,         # 预测区间粗细
                    angle = 0,
                    legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p11
coef_rotation <- coef(overall_model_ap_Stress)
vcov_rotation <- vcov(overall_model_ap_Stress)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# stressno stress    stressstress 
# "a"             "b" 














### 8.10 Native
ap_Native <- subset(ap, native %in% c("native", "no native"))
#
ap_Native$native <- droplevels(factor(ap_Native$native))
# The number of Observations and StudyID
group_summary <-ap_Native %>%
  group_by(native) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 2 × 3
# native    Observations Unique_StudyID
# <fct>            <int>          <int>
#   1 native              20              6
# 2 no native          150             40

overall_model_ap_Native <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + native, random = ~ 1 | study.id, data = ap_Native, method = "REML")
# QM and p value
summary(overall_model_ap_Native)
# Multivariate Meta-Analysis Model (k = 170; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016875.8541   4033751.7082   4033757.7082   4033767.0800   4033757.8545   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0371  0.1926     45     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 168) = 4612627.6255, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 2.6589, p-val = 0.2646
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# nativenative      -0.0114  0.0347  -0.3286  0.7424  -0.0795  0.0567    
# nativeno native    0.0177  0.0305   0.5792  0.5625  -0.0422  0.0776    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p12 <- orchard_plot(overall_model_ap_Native,
                    mod = "native",          # 截距模型用 "1"
                    group = "study.id",
                    xlab = "Effect size(lnRR)",
                    trunk.size = 0.5,        # 组均值的置信区间粗细
                    branch.size = 2,         # 预测区间粗细
                    angle = 0,
                    legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p12
coef_rotation <- coef(overall_model_ap_Native)
vcov_rotation <- vcov(overall_model_ap_Native)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# nativeno native    nativenative 
# "a"             "a" 











### 8.10 Plant.type
ap_Plant.type <- subset(ap, plant.type %in% c("crop", "herbaceous plants", "lianas", "shrub", "tree"))
#
ap_Plant.type$plant.type <- droplevels(factor(ap_Plant.type$plant.type))
# The number of Observations and StudyID
group_summary <-ap_Plant.type %>%
  group_by(plant.type) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 4 × 3
# plant.type        Observations Unique_StudyID
# <fct>                    <int>          <int>
#   1 crop                        42             14
# 2 herbaceous plants           83             22
# 3 shrub                       19              5
# 4 tree                        43             10

overall_model_ap_Plant.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + plant.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Plant.type)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016880.3216   4033760.6433   4033770.6433   4033786.6907   4033770.9823   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0296  0.1720     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 183) = 4075573.0046, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 2.7726, p-val = 0.5966
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# plant.typecrop                -0.0234  0.0439  -0.5336  0.5936  -0.1096  0.0627    
# plant.typeherbaceous plants    0.0324  0.0368   0.8805  0.3786  -0.0397  0.1045    
# plant.typeshrub                0.1190  0.0969   1.2278  0.2195  -0.0710  0.3090    
# plant.typetree                 0.0110  0.0597   0.1836  0.8543  -0.1060  0.1279    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p13 <- orchard_plot(overall_model_ap_Plant.type,
                    mod = "plant.type",          # 截距模型用 "1"
                    group = "study.id",
                    xlab = "Effect size(lnRR)",
                    trunk.size = 0.5,        # 组均值的置信区间粗细
                    branch.size = 2,         # 预测区间粗细
                    angle = 0,
                    legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p13
coef_rotation <- coef(overall_model_ap_Plant.type)
vcov_rotation <- vcov(overall_model_ap_Plant.type)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# plant.typeshrub plant.typeherbaceous plants              plant.typetree              plant.typecrop 
# "a"                         "a"                         "a"                         "a" 










### 8.11 Experimental.periods
ap_Experimental.periods <- subset(ap, experimental.periods %in% c("Short-term", "Medium-term", "Long-term"))
#
ap_Experimental.periods$experimental.periods <- droplevels(factor(ap_Experimental.periods$experimental.periods))
# The number of Observations and StudyID
group_summary <-ap_Experimental.periods %>%
  group_by(experimental.periods) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 3 × 3
# experimental.periods Observations Unique_StudyID
# <fct>                       <int>          <int>
#   1 Long-term                      20              4
# 2 Medium-term                    44             14
# 3 Short-term                    107             26

overall_model_Experimental.periods <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + experimental.periods, random = ~ 1 | study.id, data = ap_Experimental.periods, method = "REML")
# QM and p value
summary(overall_model_Experimental.periods)
# Multivariate Meta-Analysis Model (k = 171; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2016897.7564   4033795.5129   4033803.5129   4033816.0087   4033803.7583   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0338  0.1838     44     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 168) = 4060551.3881, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
#   QM(df = 3) = 0.5455, p-val = 0.9088
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    
# experimental.periodsLong-term     -0.0162  0.0998  -0.1622  0.8712  -0.2118  0.1794    
# experimental.periodsMedium-term    0.0296  0.0538   0.5493  0.5828  -0.0759  0.1350    
# experimental.periodsShort-term     0.0174  0.0373   0.4663  0.6410  -0.0557  0.0905    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p14 <- orchard_plot(overall_model_Experimental.periods,
                    mod = "experimental.periods",          # 截距模型用 "1"
                    group = "study.id",
                    xlab = "Effect size(lnRR)",
                    trunk.size = 0.5,        # 组均值的置信区间粗细
                    branch.size = 2,         # 预测区间粗细
                    angle = 0,
                    legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p14
coef_rotation <- coef(overall_model_Experimental.periods)
vcov_rotation <- vcov(overall_model_Experimental.periods)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(p)
}

# 比较
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)]

# 现在，我们按照效应大小（系数值）对组进行排序（从大到小）
effect_sizes <- coef_rotation
# 按效应大小从大到小排序
sorted_indices <- order(effect_sizes, decreasing = TRUE)
sorted_group_names <- group_names[sorted_indices]

# 重新排列p值矩阵，按照效应大小顺序
p_matrix_sorted <- p_matrix[sorted_group_names, sorted_group_names]

# 现在，使用排序后的p值矩阵计算字母
significance_letters_sorted <- multcompLetters(p_matrix_sorted)$Letters

# 输出
print(significance_letters_sorted)
# experimental.periodsMedium-term  experimental.periodsShort-term   experimental.periodsLong-term 
# "a"                             "a"                             "a" 













### 8.12 Fertilizer
ap_Fertilizer <- subset(ap, Fertilizer %in% c("Fertilizer", "No Fertilizer"))
#
ap_Fertilizer$Fertilizer <- droplevels(factor(ap_Fertilizer$Fertilizer))
# The number of Observations and StudyID
group_summary <-ap_Fertilizer %>%
  group_by(Fertilizer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)

overall_model_ap_Fertilizer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + fertilizer, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Fertilizer)
# Multivariate Meta-Analysis Model (k = 187; method: REML)
# 
# logLik       Deviance            AIC            BIC           AICc   
# -2015945.3966   4031890.7933   4031896.7933   4031906.4543   4031896.9258   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0396  0.1990     49     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 185) = 4191664.7452, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1878.8589, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     
# fertilizerfertilizer       0.0779  0.0301   2.5905  0.0096   0.0190  0.1368  ** 
#   fertilizerno fertilizer   -0.0009  0.0300  -0.0306  0.9756  -0.0598  0.0579     
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p15 <- orchard_plot(overall_model_ap_Fertilizer,
                    mod = "fertilizer",          # 截距模型用 "1"
                    group = "study.id",
                    xlab = "Effect size(lnRR)",
                    trunk.size = 0.5,        # 组均值的置信区间粗细
                    branch.size = 2,         # 预测区间粗细
                    angle = 0,
                    legend.pos = "none"      # 关闭右上角 Precision(1/SE) 图例
) +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    panel.background = element_blank(),   # 透明背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6) # 保留边框（可选）
  )
p15
coef_rotation <- coef(overall_model_ap_Fertilizer)
vcov_rotation <- vcov(overall_model_ap_Fertilizer)
# 定义多重比较函数
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2])
  z <- diff / se_diff
  2 * (1 - pnorm(abs(z)))
}

group_names <- names(coef_rotation)

# p matrix (raw)
p_matrix <- matrix(1, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }
}
diag(p_matrix) <- 1

# ---- p-adjust on upper triangle (non-NA only) ----
idx <- which(upper.tri(p_matrix) & !is.na(p_matrix), arr.ind = TRUE)
p_raw <- p_matrix[idx]
p_adj <- p.adjust(p_raw, method = "holm")   # 或 "BH"/"bonferroni"

p_matrix_adj <- p_matrix
for (k in seq_len(nrow(idx))) {
  i <- idx[k, 1]; j <- idx[k, 2]
  p_matrix_adj[i, j] <- p_adj[k]
  p_matrix_adj[j, i] <- p_adj[k]
}
diag(p_matrix_adj) <- 1

# ---- sort by effect size (safe) ----
rn <- rownames(p_matrix_adj)
sorted_group_names <- rn[order(coef_rotation[rn], decreasing = TRUE, na.last = TRUE)]
p_matrix_adj_sorted <- p_matrix_adj[sorted_group_names, sorted_group_names, drop = FALSE]

# letters
letters <- multcompLetters(p_matrix_adj_sorted)$Letters
print(letters)
# fertilizerfertilizer fertilizerno fertilizer 
# "a"                     "b" 
#patchwork
p1+p2+p3+p4+p5+p6+p7+p9+p10+p11+p12+p13+p14+p15+ plot_layout(ncol = 3)
#图片导出为23×20






library(lme4)
ap$Wr <- 1 / ap$Vi
###background character
##time
time <- ggplot(ap, aes(y=RR, x=experimental.duration)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="Experimental duration" , y="RR of phasphotase")
time
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
m1 <- lmer(RR ~ dur_z + inoculant.type + (1|study.id), weights=Wr, data=ap)
m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + (1|study.id), weights=Wr, data=ap)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# dur_z           11184   11184     1  5.3185  0.4526 0.5292
# I(dur_z^2)      14082   14082     1  6.8266  0.5699 0.4755
# inoculant.type  12747    4249     3 27.6396  0.1720 0.9144
anova(m1, m2)
#pH
sum(!is.na(ap$soil.ph))#134
ap$pH_z <- as.numeric(scale(ap$soil.ph))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~  dur_z + I(dur_z^2) + inoculant.type + pH_z + (1|study.id),
           weights = Wr, data = ap)

m2 <- lmer(RR ~  dur_z + I(dur_z^2) + inoculant.type + pH_z + I(pH_z^2) + (1|study.id),
           weights = Wr, data = ap)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            1208    1208     1 19.161  1.6715  0.211426    
# I(dur_z^2)       1282    1282     1 27.050  1.7737  0.194035    
# inoculant.type   1560     520     3 72.678  0.7196  0.543489    
# pH_z             7726    7726     1 37.989 10.6889  0.002294 ** 
# I(pH_z^2)       49254   49254     1 38.489 68.1472 4.795e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m1, m2) 
library(ggeffects)
library(ggplot2)

# 1) 用模型里的变量名做预测
pred_z <- ggpredict(m2, terms = "pH_z [all]")  # x 是 pH_z

# 2) 把 x 从 z 值换回原始 soil.ph
mu  <- attr(scale(ap$soil.ph), "scaled:center")
sdv <- attr(scale(ap$soil.ph), "scaled:scale")
pred_z$soil_pH <- pred_z$x * sdv + mu

# 3) 拐点（原始 pH）
b1 <- fixef(m2)["pH_z"]
b2 <- fixef(m2)["I(pH_z^2)"]
z_star  <- -b1/(2*b2)
pH_star <- z_star * sdv + mu
pH_star
# pH_z 
# 6.324806 
# 4) 用“模型实际使用的数据”画散点（避免 NA 行混进来）
ap_used <- model.frame(m2)                # m2 实际使用的数据行
# 注意：ap_used 里没有 soil.ph，需要从原数据取同样行号时最稳：
idx <- as.integer(rownames(ap_used))
ap_used$soil.ph <- ap$soil.ph[idx]

# 5) 作图：散点 + LMM 预测曲线 + 95%CI + 拐点
p16 <- ggplot() +
  geom_point(data = ap_used, aes(x = soil.ph, y = RR),
             color = "gray", size = 10, shape = 21) +
  geom_ribbon(data = pred_z, aes(x = soil_pH, ymin = conf.low, ymax = conf.high),
              alpha = 0.20) +
  geom_line(data = pred_z, aes(x = soil_pH, y = predicted),
            linewidth = 2) +
  geom_vline(xintercept = pH_star, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Soil pH", y = "Effect size (RR of Soil-available P)")
p16









##TC

p17 <- ggplot(ap, aes(y=RR, x=soil.tc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil TC" , y="RR of Soil-available P")
p17

ap$soil.tc_z <- as.numeric(scale(ap$soil.tc))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(
  RR ~ dur_z + I(dur_z^2) + inoculant.type + soil.tc_z + (1 | study.id),
  weights = Wr, data = ap
)

m2 <- lmer(
  RR ~ dur_z + I(dur_z^2) + inoculant.type + soil.tc_z + I(soil.tc_z^2) + (1 | study.id),
  weights = Wr, data = ap
)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF    DenDF F value Pr(>F)
# dur_z          0.61455 0.61455     1 0.073275  0.0304 0.9579
# I(dur_z^2)     0.61460 0.61460     1 0.073275  0.0304 0.9579
# inoculant.type 0.61445 0.61445     1 0.073274  0.0304 0.9579
# soil.tc_z      0.61340 0.61340     1 0.073275  0.0303 0.9579













##SOC
p18 <- ggplot(ap, aes(y=RR, x=soil.soc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil SOC" , y="RR of Soil-available P")
p18
ap$soc_z  <- as.numeric(scale(ap$soil.soc))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + I(soc_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)   # 若不显著，停止
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# dur_z          28.802  28.802     1    75  0.0668 0.7968
# I(dur_z^2)     22.190  22.190     1    75  0.0514 0.8212
# inoculant.type 72.480  36.240     2    75  0.0840 0.9195
# soc_z          34.698  34.698     1    75  0.0804 0.7775




##TN
p19 <- ggplot(ap, aes(y=RR, x=tn)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil TN" , y="RR of Soil-available P")
p19


ap$tn_z  <- as.numeric(scale(ap$tn))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tn_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tn_z + I(tn_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# dur_z            1.977   1.977     1    90  0.0054 0.9413
# I(dur_z^2)       1.406   1.406     1    90  0.0039 0.9505
# inoculant.type   1.181   0.591     2    90  0.0016 0.9984
# tn_z           127.892 127.892     1    90  0.3522 0.5544





##TP

tbl.TP<-lmer(RR~ scale(experimental.duration)+inoculant.type +scale(tp) +(1|study.id), weights=Wr, data=ap)
summary(tbl.TP)
# Number of obs: 62, groups:  study.id, 19
anova(tbl.TP)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(experimental.duration)  0.1839  0.1839     1 10.620  0.0211 0.8873
# inoculant.type               10.6877  5.3438     2 11.015  0.6119 0.5598
# scale(tp)                     0.1050  0.1050     1 44.569  0.0120 0.9132
p20 <- ggplot(ap, aes(y=RR, x=tp)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil TP" , y="RR of Soil-available P")
p20
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$tp_z  <- as.numeric(scale(ap$tp))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tp_z + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tp_z + I(tp_z^2) + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# dur_z           1.1778  1.1778     1  9.689  0.1346 0.7216
# I(dur_z^2)      2.7829  2.7829     1 11.036  0.3180 0.5841
# inoculant.type 15.5279  7.7640     2 13.197  0.8871 0.4350
# tp_z            0.3612  0.3612     1 60.805  0.0413 0.8397



##TK

ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$tk_z  <- as.numeric(scale(ap$tk))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + I(tk_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# dur_z          0.1431 0.14312     1 4.4937  0.0435 0.8438
# I(dur_z^2)     0.9623 0.96227     1 5.2500  0.2928 0.6106
# inoculant.type 3.3993 1.69967     2 5.2362  0.5171 0.6238
# tk_z           0.0143 0.01432     1 6.1519  0.0044 0.9495


##AP
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ap_z  <- as.numeric(scale(ap$ap))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ap_z + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ap_z + I(ap_z^2) + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

anova(m1, m2)

anova(m1)








##AK
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ak_z  <- as.numeric(scale(ap$ak))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ak_z + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ak_z + I(ak_z^2) + (1|study.id),
           weights=Wr, data=ap, REML=FALSE)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# dur_z          0.7162  0.7162     1  9.809  0.3930 0.5450
# I(dur_z^2)     0.0227  0.0227     1 17.420  0.0125 0.9124
# inoculant.type 6.4382  3.2191     2 13.272  1.7666 0.2088
# ak_z           1.7314  1.7314     1 15.675  0.9502 0.3445






#图片输出为10×12
p16+p17+p18+p19+p20+p22+p24+p25+plot_layout(ncol = 3)





# Unit: Bacteria_cells_kg_soil
Bacteria_cells_kg_soil <-read.csv(file.choose())
Bacteria_cells_kg_soil$Wr <- 1 / Bacteria_cells_kg_soil$Vi
mdBacteria_cells_kg_soil<-lmer(RR~ scale(final.inoculation.amount)+scale(experimental.duration) + (1|study.id), weights=Wr, data=Bacteria_cells_kg_soil)
anova(Bacteria_cells_kg_soil) 
r.squaredGLMM(mdFinalinoculationquantification)
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)
r.squaredGLMM(mdLog10)








# Unit: Bacteria_cfu_g_seed
Bacteria_cfu_g_seed <-read.csv(file.choose())
Bacteria_cfu_g_seed$Wr <- 1 / Bacteria_cfu_g_seed$Vi
##Finalinoculationquantification
mdBacteria_cfu_g_seed<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_g_seed)
anova(mdBacteria_cfu_ha_soil)
p25 <- ggplot(Bacteria_cfu_ha_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="Cfu/ha soil" , y="RR of Soil-available P")
p25
r.squaredGLMM(mdFinalinoculationquantification)
mdlBacteria_cfu_ha_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_ha_soil)
anova(mdlBacteria_cfu_ha_soil)
p26 <- ggplot(Bacteria_cfu_ha_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="log10(Cfu/ha soil)" , y="RR of Soil-available P")
p26







# Unit: Bacteria_cfu_kg_soil
Bacteria_cfu_kg_soil <-read.csv(file.choose())
Bacteria_cfu_kg_soil$Wr <- 1 / Bacteria_cfu_kg_soil$Vi
##Finalinoculationquantification
mdBacteria_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_kg_soil)
anova(mdBacteria_cfu_kg_soil)
p27 <- ggplot(Bacteria_cfu_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="Cfu/kg soil" , y="RR of Soil-available P")
p27
mdlBacteria_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_kg_soil)
anova(mdlBacteria_cfu_kg_soil)
p28 <- ggplot(Bacteria_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="log10(Cfu/kg soil)" , y="RR of Soil-available P")
p28





# Unit: Bacteria_cfu_root
Bacteria_cfu_root <-read.csv(file.choose())
Bacteria_cfu_root$Wr <- 1 / Bacteria_cfu_root$Vi
##Finalinoculationquantification
mdBacteria_cfu_root<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_root)
anova(mdBacteria_cfu_root)
r.squaredGLMM(mdFinalinoculationquantification)
mdlBacteria_cfu_root<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_root)
anova(mdlBacteria_cfu_root)






# Unit: Bacteria_cfu_seed
Bacteria_cfu_seed <-read.csv(file.choose())
Bacteria_cfu_seed$Wr <- 1 / Bacteria_cfu_seed$Vi
##Finalinoculationquantification
mdBacteria_cfu_seed<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_seed)
anova(mdBacteria_cfu_seed)
r.squaredGLMM(mdFinalinoculationquantification)
mdlBacteria_cfu_seed<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_seed)
anova(mdlBacteria_cfu_seed)










# Unit: Bacteria_spores_kg_soil
Bacteria_spores_kg_soil <-read.csv(file.choose())
Bacteria_spores_kg_soil$Wr <- 1 / Bacteria_spores_kg_soil$Vi
##Finalinoculationquantification
mdBacteria_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_spores_kg_soil)
anova(mdBacteria_spores_kg_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlBacteria_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_spores_kg_soil)
anova(mdlBacteria_spores_kg_soil)









# Unit: fungi_cfu_kg_soil
Fungi_cfu_kg_soil <-read.csv(file.choose())
Fungi_cfu_kg_soil$Wr <- 1 / Fungi_cfu_kg_soil$Vi
##Finalinoculationquantification
mdFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdFungi_cfu_kg_soil)
p29 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="cfu/kg soil" , y="RR of Soil-available P")
p29
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)
p30 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="log10(cfu/kg soil)" , y="RR of Soil-available P")
p30







# Unit: Fungi_pieces_root
Fungi_pieces_root <-read.csv(file.choose())
Fungi_pieces_root$Wr <- 1 / Fungi_pieces_root$Vi
##Finalinoculationquantification
mdFungi_pieces_root<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_pieces_root)
anova(mdFungi_cfu_kg_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)








# Unit: Fungi_spores_l_soil
Fungi_spores_l_soil <-read.csv(file.choose())
Fungi_spores_l_soil$Wr <- 1 / Fungi_spores_l_soil$Vi
##Finalinoculationquantification
mdFungi_spores_l_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_spores_l_soil)
anova(mdFungi_spores_l_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)







# Unit: Fungi_spores_seed
Fungi_spores_seed <-read.csv(file.choose())
Fungi_spores_seed$Wr <- 1 / Fungi_spores_seed$Vi
##Finalinoculationquantification
mdFungi_spores_seed<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_spores_seed)
anova(mdFungi_spores_l_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)







# Unit: Fungi_spores_seeding
Fungi_spores_seeding <-read.csv(file.choose())
Fungi_spores_seeding$Wr <- 1 / Fungi_spores_seeding$Vi
##Finalinoculationquantification
mdFungi_spores_seeding<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_spores_seeding)
anova(mdFungi_spores_l_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)







# Unit: AMF_cfu_kg_soil
AMF_cfu_kg_soil <-read.csv(file.choose())
AMF_cfu_kg_soil$Wr <- 1 / AMF_cfu_kg_soil$Vi
##Finalinoculationquantification
mdAMF_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_cfu_kg_soil)
anova(mdAMF_cfu_kg_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_cfu_kg_soil)
anova(mdlAMF_cfu_kg_soil)








# Unit: AMF_propagules_ha_soil
AMF_propagules_ha_soil <-read.csv(file.choose())
AMF_propagules_ha_soil$Wr <- 1 / AMF_propagules_ha_soil$Vi
##Finalinoculationquantification
mdAMF_propagules_ha_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_propagules_ha_soil)
anova(mdAMF_propagules_ha_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_propagules_ha_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_propagules_ha_soil)
anova(mdlAMF_propagules_ha_soil)







# Unit: AMF_propagules_kg_soil
AMF_propagules_kg_soil <-read.csv(file.choose())
AMF_propagules_kg_soil$Wr <- 1 / AMF_propagules_kg_soil$Vi
##Finalinoculationquantification
mdAMF_propagules_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_propagules_kg_soil)
anova(mdAMF_propagules_kg_soil)
p31 <- ggplot(AMF_propagules_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="propagules/kg soil" , y="RR of Soil-available P")
p31

r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_propagules_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_propagules_kg_soil)
anova(mdlAMF_propagules_kg_soil)
p32 <- ggplot(AMF_propagules_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="log10(propagules/kg soil)" , y="RR of Soil-available P")
p32









# Unit: AMF_spores_kg_soil
AMF_spores_kg_soil <-read.csv(file.choose())
AMF_spores_kg_soil$Wr <- 1 / AMF_spores_kg_soil$Vi
##Finalinoculationquantification
mdAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdAMF_spores_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration)    0.05953 0.05953     1    10  0.0218 0.8857
# scale(final.inoculation.amount) 0.93582 0.93582     1    10  0.3421 0.5716
p33 <- ggplot(AMF_spores_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="spores/kg soil" , y="RR of Soil-available P")
p33
mdlAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdlAMF_spores_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration) 0.0777  0.0777     1    10  0.0292 0.8677
# scale(log10)                 1.7149  1.7149     1    10  0.6452 0.4405
p34 <- ggplot(AMF_spores_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="log10(spores/kg soil)" , y="RR of Soil-available P")
p34



# Unit: AMF_cfu_kg_soil
AMF_cfu_kg_soil <-read.csv(file.choose())
AMF_cfu_kg_soil$Wr <- 1 /AMF_cfu_kg_soil$Vi
##Finalinoculationquantification
mdAMF_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_cfu_kg_soil)
anova(mdAMF_cfu_kg_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdlAMF_spores_kg_soil)









# Unit: AMF_spores_l_soil
AMF_spores_l_soil <-read.csv(file.choose())
AMF_spores_l_soil$Wr <- 1 /AMF_spores_l_soil$Vi
##Finalinoculationquantification
mdAMF_spores_l_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_spores_l_soil)
anova(mdAMF_spores_l_soil)
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdlAMF_spores_kg_soil)






# Unit: AMF_spores_seeding
AMF_spores_seeding <-read.csv(file.choose())
AMF_spores_seeding$Wr <- 1 / AMF_spores_seeding$Vi
##Finalinoculationquantification
mdAMF_spores_seeding<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_spores_seeding)
anova(mdAMF_spores_seeding)
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdlAMF_spores_kg_soil)

## RRpH
## RRpH
library(ggplot2)
library(ggpmisc)
library(dplyr)

# 读入三个数据
acid_soil    <- read.csv(file.choose())
N_soil       <- read.csv(file.choose())
alkline_soil <- read.csv(file.choose())

# 统一计算 Spearman（raw p、rho、n）
get_spear <- function(dat) {
  df <- dat %>% filter(!is.na(RR), !is.na(RRpH))
  ct <- suppressWarnings(cor.test(df$RR, df$RRpH, method="spearman", exact=FALSE))
  list(n = nrow(df), rho = unname(ct$estimate), p = ct$p.value)
}

s_acid <- get_spear(acid_soil)
s_N    <- get_spear(N_soil)
s_alk  <- get_spear(alkline_soil)

# Holm 矫正（✅ 与你前面一致）
p_raw_all <- c(acid = s_acid$p, N = s_N$p, alkaline = s_alk$p)
p_adj_all <- p.adjust(p_raw_all, method = "holm")

# 画图函数：保留你原来的层，只把 stat_cor 换成 annotate(矫正p)
plot_RR_RRpH <- function(dat, tag, out_pdf) {
  
  st_n   <- if (tag=="acid") s_acid$n else if (tag=="N") s_N$n else s_alk$n
  st_rho <- if (tag=="acid") s_acid$rho else if (tag=="N") s_N$rho else s_alk$rho
  
  lab_cor <- sprintf("Spearman \u03c1 = %.3f\nHolm-adjusted p = %.3g\nn = %d",
                     st_rho, p_adj_all[[tag]], st_n)
  
  p <- ggplot(dat, aes(y=RR, x=RRpH)) +
    geom_point(color="gray", size=10, shape=21) +
    geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") +
    theme_bw() +
    theme(text = element_text(size=20),
          panel.grid = element_blank()) +
    geom_hline(yintercept=0, colour="black", linewidth=0.5, linetype="dashed") +
    stat_poly_eq(
      aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
      formula = y ~ x, parse = TRUE, color="black",
      size = 5, label.x = 0.05, label.y = 0.85
    ) +
    annotate("text",
             x = -Inf, y = Inf, hjust = -0.05, vjust = 1.10,
             label = lab_cor, size = 5) +
    labs(x="RRpH", y="RR625")
  
  print(p)
  
  pdf(out_pdf, width=8, height=8)
  print(p)
  dev.off()
  
  invisible(p)
}

# 输出三张图（文件名按你原来那样）
plot_RR_RRpH(acid_soil,    tag = "acid",     out_pdf = "RRpH_acid_Holm.pdf")
plot_RR_RRpH(N_soil,       tag = "N",        out_pdf = "RRpH_N_Holm.pdf")
plot_RR_RRpH(alkline_soil, tag = "alkaline", out_pdf = "RRpH_alkline_Holm.pdf")

# packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpmisc)
library(patchwork)   # 用于拼图；没有就 install.packages("patchwork")

# -------------------------
# 0) 变量列表（去重）
# -------------------------
x_vars <- unique(c(
  "RRyield",
  "RRbiomass",
  "RRshoot_biomass",
  "RRroot_biomass",
  "RRroot_biomass",
  "RR_pleaf",
  "RR_pgrain",
  "RR_pshoot",
  "RR_proot",
  "RR_pplant"# 你重复写了，unique() 会自动去掉
))

# 检查变量是否都在数据里
miss <- setdiff(x_vars, names(ap))
if (length(miss) > 0) stop("ap 里缺少这些列：", paste(miss, collapse = ", "))

# -------------------------
# 1) 一次性计算 Spearman + Holm 矫正
# -------------------------
adjust_method <- "holm"  # ✅ 和你前面 z-test 保持一致

cor_tbl <- map_dfr(x_vars, function(x) {
  df <- ap %>% filter(!is.na(RR), !is.na(.data[[x]]))
  ct <- suppressWarnings(cor.test(df$RR, df[[x]], method = "spearman", exact = FALSE))
  tibble(
    var   = x,
    n     = nrow(df),
    rho   = unname(ct$estimate),
    p_raw = ct$p.value
  )
}) %>%
  mutate(p_adj = p.adjust(p_raw, method = adjust_method)) %>%
  arrange(p_adj)

print(cor_tbl)
write.csv(cor_tbl, "Spearman_RR_vs_traits_Holm.csv", row.names = FALSE)

# -------------------------
# 2) 生成单张图：保留你原来的风格，只把 stat_cor 换成 Holm-adjusted p 的注释
# -------------------------
make_plot <- function(x, stats_tbl = cor_tbl) {
  
  st <- stats_tbl %>% filter(var == x)
  if (nrow(st) != 1) stop("stats_tbl 里找不到该变量：", x)
  
  # 你原来那句 stat_cor 显示的内容，用 Holm-adjusted p 替换
  lab_cor <- sprintf("Spearman \u03c1 = %.3f\nHolm-adjusted p = %.3g\nn = %d",
                     st$rho, st$p_adj, st$n)
  
  ggplot(ap, aes(y = RR, x = .data[[x]])) +
    geom_point(color = "gray", size = 10, shape = 21) +
    geom_smooth(method = lm, color = "black", linewidth = 2.0,
                se = TRUE, level = 0.95, linetype = "solid") +
    theme_bw() +
    theme(
      text = element_text(family = "serif", size = 20),
      panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5, linetype = "dashed") +
    stat_poly_eq(
      aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
      formula = y ~ x, parse = TRUE, color = "black",
      size = 5,
      label.x = 0.05,
      label.y = 0.85
    ) +
    annotate("text",
             x = -Inf, y = Inf, hjust = -0.05, vjust = 1.10,
             label = lab_cor, size = 5, family = "serif") +
    labs(x = x, y = "RR")
}

plots <- map(x_vars, make_plot)

# -------------------------
# 3) 拼成多面板（例如 3 列），输出一个 PDF
# -------------------------
p_all <- wrap_plots(plots, ncol = 3)
p_all
ggsave("RR_vs_traits_HolmAdjusted.pdf", p_all, width = 16, height = 10)

library(rfPermute)
library(randomForest)
# Run the random forest and calculate the permutation test of variable importance
rf_model_perm3 <- rfPermute(RR ~ experimental.duration + soil.ph + soil.tc + soil.soc + tn + tp +
                              tk + an + ap + ak, 
                            data = ap, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
# %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# experimental.duration 24.508305   0.00990099    1.96410441         0.03960396
# soil.ph               19.736048   0.00990099    1.90497030         0.00990099
# soil.tc               13.233623   0.00990099    0.99492728         0.00990099
# tn                    12.611922   0.00990099    0.69618954         0.11881188
# soil.soc              12.259412   0.00990099    0.54044134         0.48514851
# tp                    10.146593   0.02970297    0.47586570         0.35643564
# ap                     8.653985   0.02970297    0.40515345         0.56435644
# tk                     7.877038   0.00990099    0.70072772         0.01980198
# ak                     6.008281   0.05940594    0.33951065         0.22772277
# an                     3.914552   0.09900990    0.02969814         0.73267327
# install.packages(c("ggplot2", "dplyr", "forcats"))  # 如未安装
library(ggplot2)

# 提取 importance
imp <- importance(rf_model_perm3)
imp_df <- as.data.frame(imp)
imp_df$Variable <- rownames(imp_df)

# 按 %IncMSE 从大到小排序
imp_df <- imp_df[order(imp_df$`%IncMSE`, decreasing = TRUE), ]

# 显著性分组
imp_df$Sig <- ifelse(imp_df$`%IncMSE.pval` < 0.05, "sig", "ns")

# P值标签（保留3位；如果你想显示更多位数可改 digits）
imp_df$P_lab <- paste0("P = ", formatC(imp_df$`%IncMSE.pval`, format = "f", digits = 3))

# 因子顺序（保持从左到右按降序）
imp_df$Variable <- factor(imp_df$Variable, levels = imp_df$Variable)

# 绘图（竖向柱图，更接近你示例）
library(ggplot2)

# 你前面已经得到 imp_df 的话，可以直接接着用
# imp_df 包含: Variable, `%IncMSE`, `%IncMSE.pval`, Sig, P_lab

p <- ggplot(imp_df, aes(x = Variable, y = `%IncMSE`, fill = Sig)) +
  # 减小柱宽（原来接近1，这里改成0.65）
  geom_col(width = 0.5, color = NA) +
  
  # 柱内P值（竖排）
  geom_text(
    aes(
      y = pmax(`%IncMSE` * 0.08, 0.8),
      label = P_lab
    ),
    angle = 90,
    vjust = 0.5,
    hjust = 0,
    size = 3.2,
    color = "black"
  ) +
  
  scale_fill_manual(values = c(sig = "#4A90E2", ns = "grey80")) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.2))) +
  
  # ===== 自定义Y轴坐标 =====
scale_y_continuous(
  name = "Increase in MSE (%)",
  breaks = seq(0, 60, 10),
  limits = c(0, 30),
  expand = expansion(mult = c(0, 0.05))  # 下方0，上方留5%
)+
  
  labs(x = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

p
#8*6
######################  Trials sorted by effect size
library(ggplot2)
library(dplyr)
BacterialShannon <- read.csv("ap_effects.csv", fileEncoding = "latin1")
# 95% CI + significant
df_plot <- BacterialShannon %>%
  filter(!is.na(RR), !is.na(Vi)) %>%
  mutate(
    SE = sqrt(Vi),
    CI_lower = RR - 1.96 * SE,
    CI_upper = RR + 1.96 * SE,
    EffectClass = case_when(
      CI_lower > 0 ~ "Positive",
      CI_upper < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(RR) %>%
  mutate(Index = row_number())

# 自定义颜色
effect_colors <- c("Negative" = "#F7AF34",  
                   "Neutral"  = "#dedede",  
                   "Positive" = "#448DCD")  

# 绘图
ggplot(df_plot, aes(x = Index, y = RR, color = EffectClass)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, alpha = 0.8) +
  scale_color_manual(values = effect_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Trials sorted by effect size",
       y = "Response Ratio (RR)",
       color = "Effect") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# 统计各类数量
table(df_plot$EffectClass)
# Negative  Neutral Positive 
# 41       87       59 
## 8*8
