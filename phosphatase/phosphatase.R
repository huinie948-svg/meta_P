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
ap <- read.csv("phophatase.csv",fileEncoding = "latin1")
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
# Total number of observations in the dataset: 1706 

# 2. The number of Study
unique_studyid_number <- length(unique(ap$study.id))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# umber of unique StudyID: 197 


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | study.id,  # StudyID is radom factor
                             data = ap, 
                             method = "REML")

# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161938.3650   323876.7300   323880.7300   323891.6127   323880.7371   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0783  0.2797    197     no  study.id 
# 
# Test for Heterogeneity:
#   Q(df = 1705) = 660363.8329, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# 0.2623  0.0209  12.5554  <.0001  0.2214  0.3033  *** 
#   
#   ---
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
# Test for Funnel Plot Asymmetry: z = 3.7392, p = 0.0002
# Limit Estimate (as sei -> 0):   b = 0.2860 (CI: 0.2635, 0.3084)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# 
# Average Effect Size:         0.3119 (with file drawer: 0.0010)
# Amount of Heterogeneity:     0.1274 (with file drawer: 0.1278)
# Observed Significance Level: <.0001 (with file drawer: 0.0500)
# Target Significance Level:   0.05
# 
# Fail-safe N: 445826


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
ap_method <- subset(ap, methods %in% c("acid pme", "alkaline pme", "pde"))
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
# # A tibble: 3 × 3
# methods      Observations Unique_StudyID
# <fct>               <int>          <int>
#   1 acid pme              722            116
# 2 alkaline pme          697             99
# 3 pde                   152             22
overall_model_ap_method <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + methods, random = ~ 1 | study.id, data = ap_method, method = "REML")
# QM and p value
summary(overall_model_ap_method)
# Multivariate Meta-Analysis Model (k = 1571; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -144276.7300   288553.4600   288561.4600   288582.8902   288561.4856   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0868  0.2947    178     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1568) = 535260.8903, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
#   QM(df = 3) = 22628.9797, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# methodsacid pme        0.1612  0.0231   6.9772  <.0001  0.1159  0.2065  *** 
#   methodsalkaline pme    0.3605  0.0231  15.5976  <.0001  0.3152  0.4058  *** 
#   methodspde             0.3010  0.0244  12.3264  <.0001  0.2531  0.3488  *** 
#   
#   ---
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
# methodsalkaline pme          methodspde     methodsacid pme 
# "a"                 "b"                 "c" 













### 8.2 Inoculation.location
ap_Inoculation.location <- subset(ap, inoculation.location %in% c("soils", "seedings+soils","soils+seeds","seeds+shoots", "seeds", "seedings+soils","seedings","roots"))
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
# # A tibble: 7 × 3
# inoculation.location Observations Unique_StudyID
# <fct>                       <int>          <int>
#   1 roots                          69              9
# 2 seedings                       57             10
# 3 seedings+soils                 14              2
# 4 seeds                         383             34
# 5 seeds+shoots                    8              1
# 6 soils                        1118            149
# 7 soils+seeds                    48              5

overall_model_ap_Inoculation.location <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculation.location, random = ~ 1 | study.id, data = ap_Inoculation.location, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculation.location)
# Multivariate Meta-Analysis Model (k = 1697; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127489.8215   254979.6431   254995.6431   255039.1029   254995.7287   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0899  0.2998    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1690) = 642810.1759, p-val < .0001
# 
# Test of Moderators (coefficients 1:7):
#   QM(df = 7) = 69051.6253, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# inoculation.locationroots             0.2787  0.1049   2.6565  0.0079   0.0731  0.4842   ** 
#   inoculation.locationseedings          0.2607  0.0399   6.5397  <.0001   0.1826  0.3389  *** 
#   inoculation.locationseedings+soils    0.1887  0.0379   4.9787  <.0001   0.1144  0.2630  *** 
#   inoculation.locationseeds             0.5106  0.0233  21.9178  <.0001   0.4650  0.5563  *** 
#   inoculation.locationseeds+shoots     -0.0604  0.3004  -0.2012  0.8405  -0.6492  0.5283      
# inoculation.locationsoils             0.2079  0.0230   9.0281  <.0001   0.1628  0.2531  *** 
#   inoculation.locationsoils+seeds       0.5675  0.0231  24.6092  <.0001   0.5223  0.6127  *** 
#   
#   ---
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
# inoculation.locationsoils+seeds          inoculation.locationseeds          inoculation.locationroots 
# "a"                                "b"                              "abc" 
# inoculation.locationseedings          inoculation.locationsoils inoculation.locationseedings+soils 
# "c"                                "c"                                "c" 
# inoculation.locationseeds+shoots 
# "abc" 









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
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -130075.8901   260151.7802   260161.7802   260188.9780   260161.8156   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0747  0.2733    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1702) = 525835.0231, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 63903.4482, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# inoculant.typeamf         0.2329  0.0205  11.3658  <.0001  0.1927  0.2730  *** 
#   inoculant.typebacteria    0.2260  0.0205  11.0354  <.0001  0.1859  0.2662  *** 
#   inoculant.typefungi       0.3777  0.0236  15.9795  <.0001  0.3314  0.4241  *** 
#   inoculant.typemix         0.5375  0.0205  26.2206  <.0001  0.4974  0.5777  *** 
#   
#   ---
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
# inoculant.typemix    inoculant.typefungi      inoculant.typeamf inoculant.typebacteria 
# "a"                    "b"                    "c"                    "d" 










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

overall_model_ap_Experimental.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + experimental.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Experimental.type)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161937.6689   323875.3379   323881.3379   323897.6601   323881.3520   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0787  0.2805    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 659574.4321, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 156.9042, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# experimental.typefield         0.2709  0.0464   5.8333  <.0001  0.1799  0.3619  *** 
#   experimental.typegreenhouse    0.2601  0.0235  11.0850  <.0001  0.2141  0.3061  *** 
#   
#   ---
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
# experimental.typefield experimental.typegreenhouse 
# "a"                         "a" 











### 8.5 Inoculum
ap_Inoculum <- subset(ap, inoculum %in% c("bacteria yes", "bacterial medium yes", "fungi yes", "amf yes", "amf inoculum yes", "mix yes"))
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
# A tibble: 6 × 3
# inoculum             Observations Unique_StudyID
# <fct>                       <int>          <int>
#   1 amf inoculum yes              141             35
# 2 amf yes                       191             45
# 3 bacteria yes                  570            116
# 4 bacterial medium yes           10              5
# 5 fungi yes                     110             27
# 6 mix yes                        88             18

overall_model_ap_Inoculum <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculum, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculum)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -129756.3570   259512.7140   259528.7140   259572.2164   259528.7993   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0756  0.2750    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1699) = 516531.3112, p-val < .0001
# 
# Test of Moderators (coefficients 1:7):
#   QM(df = 7) = 64548.3480, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# inoculumamf inoculum yes        0.2340  0.0209  11.1749  <.0001  0.1930  0.2751  *** 
#   inoculumamf yes                 0.2145  0.0207  10.3735  <.0001  0.1739  0.2550  *** 
#   inoculumbacteria yes            0.2254  0.0206  10.9143  <.0001  0.1849  0.2658  *** 
#   inoculumbacterial medium yes    0.3995  0.0268  14.8800  <.0001  0.3469  0.4521  *** 
#   inoculumfungal medium yes       0.3815  0.0932   4.0932  <.0001  0.1988  0.5641  *** 
#   inoculumfungi yes               0.3887  0.0239  16.2317  <.0001  0.3418  0.4357  *** 
#   inoculummix yes                 0.5527  0.0207  26.7315  <.0001  0.5122  0.5932  *** 
#   
#   ---
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
# inoculummix yes inoculumbacterial medium yes            inoculumfungi yes    inoculumfungal medium yes 
# "a"                          "b"                          "b"                       "abcd" 
# inoculumamf inoculum yes         inoculumbacteria yes              inoculumamf yes 
# "c"                          "c"                          "d" 










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

overall_model_ap_Inoculant.quantity <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculant.quantity, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculant.quantity)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -154747.5308   309495.0617   309501.0617   309517.3839   309501.0758   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0832  0.2884    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 620279.6032, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 14535.9050, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# inoculant.quantityno     0.3142  0.0215  14.6017  <.0001  0.2720  0.3564  *** 
#   inoculant.quantityyes    0.1286  0.0215   5.9698  <.0001  0.0864  0.1708  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ap_Inoculant.quantity)
vcov_rotation <- vcov(overall_model_ap_Inoculant.quantity)
# define
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
# inoculant.quantityno inoculant.quantityyes 
# "a"                   "b" 










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
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161934.6579   323869.3158   323881.3158   323913.9496   323881.3654   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0791  0.2813    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1701) = 652836.5051, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 158.1529, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub      
# soil.origincropland      0.2661  0.0413  6.4396  <.0001   0.1851  0.3471  *** 
#   soil.originforest        0.2105  0.0942  2.2343  0.0255   0.0258  0.3951    * 
#   soil.origingrassland     0.1446  0.1995  0.7246  0.4687  -0.2465  0.5357      
# soil.originother         0.2591  0.0262  9.8835  <.0001   0.2077  0.3105  *** 
#   soil.originplantation    0.3897  0.1061  3.6721  0.0002   0.1817  0.5978  *** 
#   
#   ---
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
# soil.originplantation   soil.origincropland      soil.originother     soil.originforest  soil.origingrassland 
# "a"                   "a"                   "a"                   "a"                   "a"













### 8.8 Sandy.or.not.sandy
ap_Sandy.or.not.sandy <- subset(ap, sandy.or.not.sandy %in% c("sandy", "not sandy"))
#
ap_Sandy.or.not.sandy$sandy.or.not.sandy <- droplevels(factor(ap_Sandy.or.not.sandy$sandy.or.not.sandy))
# The number of Observations and StudyID
group_summary <- ap_Sandy.or.not.sandy %>%
  group_by(sandy.or.not.sandy) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# # A tibble: 2 × 3
# sandy.or.not.sandy Observations Unique_StudyID
# <fct>                     <int>          <int>
#   1 not sandy                   318             30
# 2 sandy                       370             40

overall_model_ap_Sandy.or.not.sandy <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + sandy.or.not.sandy, random = ~ 1 | study.id, data = ap_Sandy.or.not.sandy, method = "REML")
# QM and p value
summary(overall_model_ap_Sandy.or.not.sandy)
# Multivariate Meta-Analysis Model (k = 688; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -31386.8873   62773.7746   62779.7746   62793.3672   62779.8098   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1000  0.3163     65     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 686) = 152245.6334, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 53.4827, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# sandy.or.not.sandynot sandy    0.2300  0.0434  5.3036  <.0001  0.1450  0.3150  *** 
#   sandy.or.not.sandysandy        0.3055  0.0424  7.1991  <.0001  0.2224  0.3887  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract model coefficients and covariance matrix
p8 <- orchard_plot(overall_model_ap_Sandy.or.not.sandy,
                   mod = "sandy.or.not.sandy",          # 截距模型用 "1"
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
p8
coef_rotation <- coef(overall_model_ap_Sandy.or.not.sandy)
vcov_rotation <- vcov(overall_model_ap_Sandy.or.not.sandy)
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
# sandy.or.not.sandysandy sandy.or.not.sandynot sandy 
# "a"                         "b" 






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
# extraction.part Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 bulk                     804             88
# 2 rhizosphere              899            115

overall_model_ap_extraction.part <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + extraction.part, random = ~ 1 | study.id, data = ap_extraction.part, method = "REML")
# QM and p value
summary(overall_model_ap_extraction.part)
# Multivariate Meta-Analysis Model (k = 1703; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161440.9949   322881.9898   322887.9898   322904.3067   322888.0039   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0801  0.2830    196     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1701) = 531976.7706, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1154.4897, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# extraction.partbulk           0.2897  0.0212  13.6655  <.0001  0.2481  0.3312  *** 
#   extraction.partrhizosphere    0.2425  0.0212  11.4470  <.0001  0.2010  0.2841  *** 
#   
#   ---
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
# extraction.partbulk extraction.partrhizosphere 
# "a"                        "b" 

library(dplyr)
library(metafor)

phos_dat <- subset(
  ap,
  extraction.part %in% c("rhizosphere", "bulk") &
    methods %in% c("acid pme", "alkaline pme", "pde")
)

phos_dat$extraction.part <- droplevels(factor(phos_dat$extraction.part))
phos_dat$methods <- droplevels(factor(phos_dat$methods))

# 查看每组样本量
group_summary <- phos_dat %>%
  group_by(methods, extraction.part) %>%
  summarise(
    Observations = n(),
    Unique_StudyID = n_distinct(study.id),
    .groups = "drop"
  )

print(group_summary)

# acid pme
dat_acid <- subset(phos_dat, methods == "acid pme")

model_acid_mean <- rma.mv(
  yi = RR,
  V = Vi,
  mods = ~ 0 + extraction.part,
  random = ~ 1 | study.id,
  data = dat_acid,
  method = "REML"
)

summary(model_acid_mean)
# Multivariate Meta-Analysis Model (k = 722; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -24771.4753   49542.9507   49548.9507   49562.6884   49548.9842   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0930  0.3049    116     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 720) = 154427.9096, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 107.7502, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# extraction.partbulk           0.1854  0.0323  5.7388  <.0001  0.1221  0.2487  *** 
#   extraction.partrhizosphere    0.3063  0.0310  9.8788  <.0001  0.2455  0.3671  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p_model_acid_mean <- orchard_plot(model_acid_mean,
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
p_model_acid_mean
# 9*6
coef_rotation <- coef(model_acid_mean)
vcov_rotation <- vcov(model_acid_mean)
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
# extraction.partrhizosphere        extraction.partbulk 
# "a"                        "b" 

# alkaline pme
dat_alk <- subset(phos_dat, methods == "alkaline pme")

model_alk_mean <- rma.mv(
  yi = RR,
  V = Vi,
  mods = ~ 0 + extraction.part,
  random = ~ 1 | study.id,
  data = dat_alk,
  method = "REML"
)

summary(model_alk_mean)
# Multivariate Meta-Analysis Model (k = 694; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -112795.9106   225591.8212   225597.8212   225611.4400   225597.8561   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0838  0.2894     98     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 692) = 311855.4438, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1055.6228, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# extraction.partbulk           0.2831  0.0307  9.2259  <.0001  0.2230  0.3433  *** 
#   extraction.partrhizosphere    0.2362  0.0307  7.6977  <.0001  0.1760  0.2963  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p_model_alk_mean <- orchard_plot(model_alk_mean,
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
p_model_alk_mean
# 9*6
coef_rotation <- coef(model_alk_mean)
vcov_rotation <- vcov(model_alk_mean)
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
# extraction.partbulk extraction.partrhizosphere 
# "a"                        "b" 
p_model_alk_mean <- orchard_plot(model_alk_mean,
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
p_model_alk_mean
coef_rotation <- coef(model_alk_mean)
vcov_rotation <- vcov(model_alk_mean)
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
# extraction.partbulk extraction.partrhizosphere 
# "a"                        "b" 










# pde
dat_pde <- subset(phos_dat, methods == "pde")

model_pde_mean <- rma.mv(
  yi = RR,
  V = Vi,
  mods = ~ 0 + extraction.part,
  random = ~ 1 | study.id,
  data = dat_pde,
  method = "REML"
)

summary(model_pde_mean)
# Multivariate Meta-Analysis Model (k = 152; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -360.1443   720.2885   726.2885   735.3204   726.4529   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1156  0.3400     22     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 150) = 4786.2986, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 244.4900, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub      
# extraction.partbulk           0.5019  0.0755  6.6440  <.0001   0.3538  0.6499  *** 
#   extraction.partrhizosphere    0.0358  0.0772  0.4634  0.6431  -0.1155  0.1871      
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p_model_pde_mean <- orchard_plot(model_pde_mean,
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
p_model_pde_mean
coef_rotation <- coef(model_pde_mean)
vcov_rotation <- vcov(model_pde_mean)
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
# extraction.partbulk extraction.partrhizosphere 
# "a"                        "b" 

















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
# # A tibble: 2 × 3
# sterilization   Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 sterilization            672             59
# 2 unsterilization         1034            142

overall_model_ap_Sterilization <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + sterilization, random = ~ 1 | study.id, data = ap_Sterilization, method = "REML")
# QM and p value
summary(overall_model_ap_Sterilization)

# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161392.1739   322784.3478   322790.3478   322806.6700   322790.3620   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0918  0.3029    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 633648.1682, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1230.8191, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# sterilizationsterilization      0.5838  0.0246  23.7776  <.0001  0.5357  0.6319  *** 
#   sterilizationunsterilization    0.1298  0.0229   5.6687  <.0001  0.0850  0.1747  *** 
#   
#   ---
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
# sterilizationsterilization sterilizationunsterilization 
# "a"                          "b" 











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


overall_model_ap_Stress <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + stress, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Stress)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161821.0325   323642.0650   323648.0650   323664.3872   323648.0791   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0786  0.2803    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 660264.7301, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 395.8552, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# stressno stress    0.2672  0.0209  12.7642  <.0001  0.2262  0.3082  *** 
#   stressstress       0.2180  0.0211  10.3166  <.0001  0.1766  0.2594  *** 
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
#   1 native              22              4
# 2 no native          122             20

overall_model_ap_Native <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + native, random = ~ 1 | study.id, data = ap_Native, method = "REML")
# QM and p value
summary(overall_model_ap_Native)
# Multivariate Meta-Analysis Model (k = 1364; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -62121.7857  124243.5715  124249.5715  124265.2216  124249.5891   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0805  0.2838    180     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1362) = 362135.1084, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 1636.3585, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# nativenative       0.3480  0.0223  15.6246  <.0001  0.3044  0.3917  *** 
#   nativeno native    0.2504  0.0222  11.2996  <.0001  0.2070  0.2939  *** 
#   
#   ---
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
# nativenative nativeno native 
# "a"             "b" 











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
# # A tibble: 5 × 3
# plant.type        Observations Unique_StudyID
# <fct>                    <int>          <int>
#   1 crop                       490             52
# 2 herbaceous plants          840            101
# 3 lianas                      15              3
# 4 shrub                      151             21
# 5 tree                       166             26

overall_model_ap_Plant.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + plant.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Plant.type)
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -161932.8073   323865.6145   323877.6145   323910.2484   323877.6641   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0776  0.2786    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1701) = 516846.3418, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 169.3465, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# plant.typecrop                 0.2603  0.0240  10.8655  <.0001  0.2133  0.3072  *** 
#   plant.typeherbaceous plants    0.2347  0.0228  10.2787  <.0001  0.1899  0.2794  *** 
#   plant.typeliaNAs               0.3680  0.1619   2.2728  0.0230  0.0507  0.6854    * 
#   plant.typeshrub                0.3211  0.0512   6.2763  <.0001  0.2208  0.4214  *** 
#   plant.typetree                 0.3244  0.0382   8.4960  <.0001  0.2496  0.3993  *** 
#   
#   ---
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
# plant.typeliaNAs              plant.typetree             plant.typeshrub              plant.typecrop plant.typeherbaceous plants 
# "a"                         "a"                        "a"                         "a"                         "a" 










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
#   1 Long-term                     122             15
# 2 Medium-term                   424             63
# 3 Short-term                    829            110

overall_model_Experimental.periods <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + experimental.periods, random = ~ 1 | study.id, data = ap_Experimental.periods, method = "REML")
# QM and p value
summary(overall_model_Experimental.periods)
# Multivariate Meta-Analysis Model (k = 1375; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -116446.0653   232892.1306   232900.1306   232921.0267   232900.1599   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1000  0.3163    173     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1372) = 351461.2698, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
#   QM(df = 3) = 78896.6431, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# experimental.periodsLong-term      0.4302  0.0268  16.0288  <.0001   0.3776  0.4828  *** 
#   experimental.periodsMedium-term    0.0474  0.0268   1.7670  0.0772  -0.0052  0.0999    . 
# experimental.periodsShort-term     0.3586  0.0257  13.9331  <.0001   0.3082  0.4091  *** 
#   
#   ---
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
# experimental.periodsLong-term  experimental.periodsShort-term experimental.periodsMedium-term 
# "a"                             "b"                             "c" 













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
# Multivariate Meta-Analysis Model (k = 1706; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -160995.8501   321991.7002   321997.7002   322014.0224   321997.7143   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.0822  0.2867    197     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1704) = 539222.6074, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 2041.0596, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# fertilizerfertilizer       0.3109  0.0214  14.5122  <.0001  0.2689  0.3529  *** 
#   fertilizerno fertilizer    0.2471  0.0214  11.5508  <.0001  0.2052  0.2891  *** 
#   
#   ---
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
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+ plot_layout(ncol = 3)
#图片导出为23×20






library(lme4)
ap$Wr <- 1 / ap$Vi
###background character
##time

ap$dur_z <- as.numeric(scale(ap$experimental.duration))
m1 <- lmer(RR ~ dur_z + inoculant.type + (1|study.id), weights=Wr, data=ap)
m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + (1|study.id), weights=Wr, data=ap)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# dur_z           31888   31888     1 1366.17 193.436 < 2.2e-16 ***
# I(dur_z^2)       5479    5479     1 1368.85  33.234 1.008e-08 ***
# inoculant.type  42189   14063     3  428.79  85.309 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(m1, m2)
# m1: RR ~ dur_z + inoculant.type + (1 | study.id)
# m2: RR ~ dur_z + I(dur_z^2) + inoculant.type + (1 | study.id)
# npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
# m1    7 3764.9 3801.5 -1875.5    3750.9                         
# m2    8 3733.9 3775.7 -1859.0    3717.9 32.986  1  9.284e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b1 <- fixef(m2)["dur_z"]
b2 <- fixef(m2)["I(dur_z^2)"]

z_star <- -b1/(2*b2)

mu  <- attr(scale(ap$experimental.duration), "scaled:center")
sdv <- attr(scale(ap$experimental.duration), "scaled:scale")

t_star <- mu + z_star*sdv
t_star
# dur_z 
# 1297.83 
library(ggeffects)
library(ggplot2)

pred <- ggpredict(m2, terms = "dur_z [all]")

mu  <- attr(scale(ap$experimental.duration), "scaled:center")
sdv <- attr(scale(ap$experimental.duration), "scaled:scale")
pred$duration_day <- pred$x * sdv + mu

# 模型实际使用的数据行（避免 NA 不一致）
ap_used <- model.frame(m2)
idx <- as.integer(rownames(ap_used))
ap_used$experimental.duration <- ap$experimental.duration[idx]

ggplot() +
  geom_point(data=ap_used, aes(x=experimental.duration, y=RR),
             color="gray", size=10,shape=21) +
  geom_ribbon(data=pred, aes(x=duration_day, ymin=conf.low, ymax=conf.high),
              alpha=0.2) +
  geom_line(data=pred, aes(x=duration_day, y=predicted), linewidth=1.2) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  labs(x="Experimental duration (days)", y="Effect size (RR of available P)")
sum(!is.na(ap$experimental.duration)) #1375











#pH

ap$pH_z <- as.numeric(scale(ap$soil.ph))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + (1|study.id),
           weights = Wr, data = ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + I(pH_z^2) + (1|study.id),
           weights = Wr, data = ap)

anova(m1, m2) 
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            7647  7646.8     1 944.35 48.8667 5.176e-12 ***
# I(dur_z^2)        513   513.4     1 941.99  3.2811    0.0704 .  
# inoculant.type  37271 12423.6     3 353.60 79.3926 < 2.2e-16 ***
# pH_z             2877  2877.3     1 597.31 18.3871 2.102e-05 ***
# I(pH_z^2)        3394  3394.3     1 525.50 21.6912 4.062e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
# 5.911301  
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
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
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
anova(m2)



##SOC
ap$soc_z  <- as.numeric(scale(ap$soil.soc))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + I(soc_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)   # 若不显著，停止
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            9524  9524.5     1 581.00 37.9290  1.37e-09 ***
# I(dur_z^2)        376   376.5     1 582.04  1.4992    0.2213    
# inoculant.type  39441 13147.0     3 180.13 52.3548 < 2.2e-16 ***
# soc_z               0     0.1     1  46.65  0.0003    0.9863    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
library(ggeffects)
library(ggplot2)

# 1) 用模型里的变量名做预测
pred_z <- ggpredict(m2, terms = "soc_z [all]")   # x 是 soc_z

# 2) 把 x 从 z 值换回原始 soil.soc
mu  <- mean(ap$soil.soc, na.rm = TRUE)
sdv <- sd(ap$soil.soc, na.rm = TRUE)
pred_z$SOC <- pred_z$x * sdv + mu

# 3) 拐点（原始 SOC）
b1 <- fixef(m2)["soc_z"]
b2 <- fixef(m2)["I(soc_z^2)"]
z_star  <- -b1 / (2 * b2)
soc_star <- z_star * sdv + mu
soc_star

# 4) 用模型实际使用的数据画散点
ap_used <- model.frame(m2)
idx <- as.integer(rownames(ap_used))
ap_used$soil.soc <- ap$soil.soc[idx]

# 5) 作图
p18 <- ggplot() +
  geom_point(data = ap_used, aes(x = soil.soc, y = RR),
             color = "gray", size = 10, shape = 21) +
  geom_ribbon(data = pred_z, aes(x = SOC, ymin = conf.low, ymax = conf.high),
              alpha = 0.20) +
  geom_line(data = pred_z, aes(x = SOC, y = predicted),
            linewidth = 1.2) +
  geom_vline(xintercept = soc_star, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Soil SOC", y = "Effect size (RR of Soil-available P)")

p18







##TN
ap$tn_z  <- as.numeric(scale(ap$tn))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tn_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tn_z + I(tn_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# dur_z          14.1443 14.1443     1 45.673  0.6412 0.4274
# I(dur_z^2)      0.4780  0.4780     1 53.636  0.0217 0.8835
# inoculant.type 16.1843  5.3948     3 90.915  0.2446 0.8650
# tn_z            0.0917  0.0917     1 32.574  0.0042 0.9490
p19 <- ggplot(ap, aes(y=RR, x=soil.soc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
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
  labs(x="Soil TN" , y="RR of Soil-available P")
p19







##TP
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$tp_z  <- as.numeric(scale(ap$tp))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tp_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tp_z + I(tp_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z               2     2.4     1 202.10  0.0078   0.92959    
# I(dur_z^2)       2910  2910.3     1 211.75  9.3477   0.00252 ** 
# inoculant.type  35541 11846.9     3  87.03 38.0509 8.318e-16 ***
# tp_z                3     3.4     1 103.35  0.0110   0.91663    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p20 <- ggplot(ap, aes(y=RR, x=soil.soc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
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
  labs(x="Soil TP" , y="RR of Soil-available P")
p20






##OP
p21 <- ggplot(ap, aes(y=RR, x=organic.p)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="organic.p" , y="RR of Soil-available P")
p21
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$op_z  <- as.numeric(scale(ap$organic.p))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + op_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + op_z + I(op_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# dur_z          0.0652  0.0652     1    24  0.0049 0.9449
# I(dur_z^2)     0.0147  0.0147     1    24  0.0011 0.9738
# inoculant.type 0.6893  0.2298     3    24  0.0172 0.9968
# op_z           9.8683  9.8683     1    24  0.7389 0.3985













##TK
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$tk_z  <- as.numeric(scale(ap$tk))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + I(tk_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# dur_z           0.0923  0.0923     1 52.899  0.0121 0.9130
# I(dur_z^2)      0.2769  0.2769     1 52.524  0.0362 0.8499
# inoculant.type  6.9765  2.3255     3  8.785  0.3036 0.8222
# tk_z           18.3419 18.3419     1  6.081  2.3949 0.1720
# I(tk_z^2)      19.5357 19.5357     1  5.695  2.5508 0.1640
# 1) 预测
pred_z <- ggpredict(m2, terms = "tk_z [all]")   # x 是 tk_z

# 2) 换回原始 tk
mu  <- mean(ap$tk, na.rm = TRUE)
sdv <- sd(ap$tk, na.rm = TRUE)
pred_z$TK <- pred_z$x * sdv + mu

# 3) 拐点
b1 <- fixef(m2)["tk_z"]
b2 <- fixef(m2)["I(tk_z^2)"]
z_star <- -b1 / (2 * b2)
tk_star <- z_star * sdv + mu
tk_star

# 4) 模型实际使用的数据
ap_used <- model.frame(m2)
idx <- as.integer(rownames(ap_used))
ap_used$tk <- ap$tk[idx]

# 5) 作图
p22 <- ggplot() +
  geom_point(data = ap_used, aes(x = tk, y = RR),
             color = "gray", size = 10, shape = 21) +
  geom_ribbon(data = pred_z, aes(x = TK, ymin = conf.low, ymax = conf.high),
              alpha = 0.20) +
  geom_line(data = pred_z, aes(x = TK, y = predicted),
            linewidth = 1.2) +
  geom_vline(xintercept = tk_star, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Soil TK", y = "Effect size (RR of Soil-available P)")

p22







##AN
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$an_z  <- as.numeric(scale(ap$an))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + an_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + an_z + I(an_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# dur_z            2.72   2.724     1 21.8748  0.0982 0.75703  
# I(dur_z^2)       8.46   8.461     1 17.6584  0.3049 0.58777  
# inoculant.type 331.03 110.342     3 14.6511  3.9758 0.02928 *
# an_z            14.05  14.048     1  4.7931  0.5062 0.50991  
# I(an_z^2)        8.15   8.152     1  5.7739  0.2937 0.60811  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p23 <- ggplot(ap, aes(y=RR, x=soil.soc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
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
  labs(x="Soil AN" , y="RR of Soil-available P")
p23




##AP
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ap_z  <- as.numeric(scale(ap$ap))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ap_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ap_z + I(ap_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            9451  9450.5     1 578.01 36.7333 2.445e-09 ***
# I(dur_z^2)        391   390.7     1 577.17  1.5186    0.2183    
# inoculant.type  39827 13275.7     3 182.99 51.6015 < 2.2e-16 ***
# ap_z               19    19.5     1 236.12  0.0757    0.7835    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p24 <- ggplot(ap, aes(y=RR, x=soil.soc)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
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
  labs(x="Soil AP" , y="RR of Soil-available P")
p24








##AK
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ak_z  <- as.numeric(scale(ap$ak))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ak_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ak_z + I(ak_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF   DenDF  F value    Pr(>F)    
# dur_z          9976.3  9976.3     1  62.015 271.2401 < 2.2e-16 ***
# I(dur_z^2)     5625.5  5625.5     1  54.029 152.9494 < 2.2e-16 ***
# inoculant.type  375.8   125.3     3 170.433   3.4058   0.01898 *  
# ak_z           1699.8  1699.8     1 174.460  46.2139 1.618e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p25 <- ggplot(ap, aes(y=RR, x=ak)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
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
  labs(x="Soil AK" , y="RR of Soil-available P")
p25

#图片输出为11×12
p16+p18+p19+p20+p21+p22+p23+p24+p25+plot_layout(ncol = 3)





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
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)    
# scale(experimental.duration)    4840.4  4840.4     1 101.857 79.5460 2.02e-14 ***
# scale(final.inoculation.amount)    6.3     6.3     1  96.695  0.1032   0.7487    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p27 <- ggplot(Bacteria_cfu_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Cfu/kg soil" , y="RR of Soil-available P")
p27
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.003183279 0.3380941
##Log10
mdlBacteria_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_kg_soil)
anova(mdlBacteria_cfu_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(experimental.duration) 4818.6  4818.6     1 101.661  79.677 1.973e-14 ***
# scale(log10)                   36.0    36.0     1  57.213   0.596    0.4433    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p28 <- ggplot(Bacteria_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
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
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(experimental.duration)    25.768  25.768     1 0.93072  0.8877 0.5281
# scale(final.inoculation.amount) 25.611  25.611     1 0.93830  0.8822 0.5280
p29 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="cfu/kg soil" , y="RR of Soil-available P")
p29
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(experimental.duration)  65.507  65.507     1 1.3822  2.2566 0.3215
# scale(log10)                 101.093 101.093     1 1.0830  3.4824 0.2983
p30 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
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
# Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration)    28.9026 28.9026     1     2  0.5584 0.5328
# scale(final.inoculation.amount)  0.0321  0.0321     1     0  0.0006 1.0000
p31 <- ggplot(AMF_propagules_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="propagules/kg soil" , y="RR of Soil-available P")
p31

r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_propagules_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_propagules_kg_soil)
anova(mdlAMF_propagules_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration) 28.9026 28.9026     1     2  0.5584 0.5328
# scale(log10)                  0.0223  0.0223     1     0  0.0004 1.0000
p32 <- ggplot(AMF_propagules_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="log10(propagules/kg soil)" , y="RR of Soil-available P")
p32









# Unit: AMF_spores_kg_soil
AMF_spores_kg_soil <-read.csv(file.choose())
AMF_spores_kg_soil$Wr <- 1 / AMF_spores_kg_soil$Vi
##Finalinoculationquantification
mdAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdAMF_spores_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(experimental.duration)    7.0802  7.0802     1 13.492  1.0050 0.3337
# scale(final.inoculation.amount) 4.4748  4.4748     1 15.957  0.6352 0.4372
p33 <- ggplot(AMF_spores_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="spores/kg soil" , y="RR of Soil-available P")
p33
r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_spores_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_spores_kg_soil)
anova(mdlAMF_spores_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(experimental.duration) 13.279  13.279     1 13.5154  1.8840 0.1922
# scale(log10)                 10.636  10.636     1  8.4643  1.5091 0.2523
p34 <- ggplot(AMF_spores_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="log10(spores/kg soil)" , y="RR of Soil-available P")
p34






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
  "RRleaf_biomass",
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
rf_model_perm3 <- rfPermute(RR ~ experimental.duration + soil.ph + soil.tc + soil.soc + tn + tp + organic.p +
                              tk + an + ap + ak, 
                            data = ap, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
# %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# soil.ph               73.055441   0.00990099    44.5521415         0.00990099
# experimental.duration 57.983056   0.00990099    30.9043527         0.00990099
# ak                    38.895126   0.00990099     7.1324831         0.00990099
# ap                    38.831565   0.00990099    12.5873430         0.00990099
# soil.soc              35.236840   0.00990099    12.5810015         0.00990099
# an                    22.877720   0.00990099     1.9877559         0.05940594
# tp                    21.916791   0.00990099     4.1441267         0.02970297
# tn                    21.233639   0.00990099     3.9843331         0.05940594
# tk                    13.783772   0.00990099     1.0114450         0.36633663
# soil.tc               11.858534   0.00990099     0.4917007         0.38613861
# organic.p             -0.439027   0.61386139     0.5403443         0.41584158

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
  limits = c(-5, 80),
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
# 91      522     1093 
## 8*8

