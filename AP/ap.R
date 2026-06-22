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
ap <- read.csv("ap.csv",fileEncoding = "latin1")
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
# Total number of observations in the dataset: 1114 

# 2. The number of Study
unique_studyid_number <- length(unique(ap$study.id))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# umber of unique StudyID: 205 


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | study.id,  # StudyID is radom factor
                             data = ap, 
                             method = "REML")

# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127445.0767   254890.1533   254894.1533   254904.1830   254894.1642   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1421  0.3770    205     no  study.id 
# 
# Test for Heterogeneity:
#   Q(df = 1113) = 5637721.6833, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# 0.2300  0.0269  8.5589  <.0001  0.1773  0.2826  *** 
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
# Test for Funnel Plot Asymmetry: z = 3.0115, p = 0.0026
# Limit Estimate (as sei -> 0):   b = 0.2498 (CI: 0.2067, 0.2928)
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# 
# Average Effect Size:         0.2931 (with file drawer: 0.0037)
# Amount of Heterogeneity:     0.2961 (with file drawer: 0.2972)
# Observed Significance Level: <.0001 (with file drawer: 0.0500)
# Target Significance Level:   0.05
# 
# Fail-safe N: 83826


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
ap_method <- subset(ap, methods %in% c("bray", "olsen", "mehlich", "resin", "water-extractable p"))
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
# # A tibble: 5 × 3
# methods             Observations Unique_StudyID
# <fct>                      <int>          <int>
#   1 bray                         135             28
# 2 mehlich                       36              3
# 3 olsen                        685            122
# 4 resin                         21              4
# 5 water-extractable p            5              1
overall_model_ap_method <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + methods, random = ~ 1 | study.id, data = ap_method, method = "REML")
# QM and p value
summary(overall_model_ap_method)
# Multivariate Meta-Analysis Model (k = 882; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -104980.1236   209960.2472   209972.2472   210000.9062   209972.3437   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1831  0.4278    156     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 877) = 508083.0877, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 44.1587, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub      
# methodsbray                   0.2970  0.0825  3.6015  0.0003   0.1354  0.4587  *** 
#   methodsmehlich                0.1335  0.2471  0.5404  0.5890  -0.3507  0.6178      
# methodsolsen                  0.2167  0.0390  5.5526  <.0001   0.1402  0.2932  *** 
#   methodsresin                  0.2128  0.0527  4.0351  <.0001   0.1094  0.3161  *** 
#   methodswater-extractable p    0.1937  0.4294  0.4512  0.6519  -0.6479  1.0354      
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
# methodsbray               methodsolsen               methodsresin methodswater-extractable p             methodsmehlich 
# "a"                        "a"                        "a"                        "a"                        "a" 













### 8.2 Inoculation.location
ap_Inoculation.location <- subset(ap, inoculation.location %in% c("roots", "roots+soils","seedings","seedings+soils", "seeds", "soils+seeds","soils"))
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
#   1 roots                          67             11
# 2 roots+soils                    25              3
# 3 seedings                       68             11
# 4 seedings+soils                 15              3
# 5 seeds                         148             31
# 6 soils                         770            153
# 7 soils+seeds                    13              4

overall_model_ap_Inoculation.location <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculation.location, random = ~ 1 | study.id, data = ap_Inoculation.location, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculation.location)
# Multivariate Meta-Analysis Model (k = 1106; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -126423.8988   252847.7976   252863.7976   252903.8148   252863.9297   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1450  0.3808    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1099) = 1456715.7825, p-val < .0001
# 
# Test of Moderators (coefficients 1:7):
#   QM(df = 7) = 1957.9658, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# inoculation.locationroots             0.2285  0.0275   8.3231  <.0001  0.1747  0.2823  *** 
#   inoculation.locationroots+soils       0.1076  0.0276   3.9004  <.0001  0.0535  0.1617  *** 
#   inoculation.locationseedings          0.1172  0.0348   3.3668  0.0008  0.0490  0.1854  *** 
#   inoculation.locationseedings+soils    0.3701  0.0337  10.9797  <.0001  0.3040  0.4361  *** 
#   inoculation.locationseeds             0.3479  0.0339  10.2657  <.0001  0.2815  0.4144  *** 
#   inoculation.locationsoils             0.2127  0.0274   7.7552  <.0001  0.1590  0.2665  *** 
#   inoculation.locationsoils+seeds       0.2134  0.0679   3.1420  0.0017  0.0803  0.3465   ** 
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
# inoculation.locationseedings+soils          inoculation.locationseeds          inoculation.locationroots 
# "abcd"                                "a"                                "b" 
# inoculation.locationsoils    inoculation.locationseeds+soils       inoculation.locationseedings 
# "c"                             "abcd"                                "d" 
# inoculation.locationroots+soils 
# "d" 








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
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127203.8347   254407.6694   254417.6694   254442.7300   254417.7237   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1421  0.3769    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1110) = 2376700.6906, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
#   QM(df = 4) = 570.8623, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# inoculant.typeamf         0.2213  0.0269  8.2273  <.0001  0.1685  0.2740  *** 
#   inoculant.typebacteria    0.2428  0.0269  9.0308  <.0001  0.1901  0.2955  *** 
#   inoculant.typefungi       0.2031  0.0280  7.2567  <.0001  0.1483  0.2580  *** 
#   inoculant.typemix         0.2037  0.0270  7.5486  <.0001  0.1508  0.2566  *** 
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
# inoculant.typebacteria      inoculant.typeamf      inoculant.typemix    inoculant.typefungi 
# "a"                    "b"                    "c"                    "bc" 










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
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127419.8667   254839.7334   254845.7334   254860.7751   254845.7551   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1532  0.3914    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1112) = 5630684.0260, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 118.8143, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# experimental.typefield         0.4138  0.0380  10.8932  <.0001  0.3394  0.4883  *** 
#   experimental.typegreenhouse    0.1699  0.0291   5.8361  <.0001  0.1129  0.2270  *** 
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
# "a"                         "b" 











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

overall_model_ap_Inoculum <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculum, random = ~ 1 | study.id, data = ap_Inoculum, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculum)
# Multivariate Meta-Analysis Model (k = 1110; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -126528.3973   253056.7945   253070.7945   253105.8414   253070.8967   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1416  0.3763    204     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1104) = 2351757.2383, p-val < .0001
# 
# Test of Moderators (coefficients 1:6):
#   QM(df = 6) = 1858.1245, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# inoculumamf inoculum yes        0.2722  0.0270  10.0953  <.0001   0.2193  0.3250  *** 
#   inoculumamf yes                 0.0341  0.0274   1.2446  0.2133  -0.0196  0.0879      
# inoculumbacteria yes            0.2892  0.0270  10.7305  <.0001   0.2364  0.3420  *** 
#   inoculumbacterial medium yes    0.4108  0.0361  11.3879  <.0001   0.3401  0.4815  *** 
#   inoculumfungi yes               0.2246  0.0281   7.9964  <.0001   0.1696  0.2797  *** 
#   inoculummix yes                 0.2265  0.0270   8.3839  <.0001   0.1736  0.2795  *** 
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
# inoculumbacterial medium yes         inoculumbacteria yes     inoculumamf inoculum yes              inoculummix yes            inoculumfungi yes 
# "a"                          "b"                          "c"                          "d"                          "d" 
# inoculumamf yes 
# "e" 










### 8.6 Inoculant.quantity
ap_Inoculant.quantity <- subset(ap, Inoculant.quantity %in% c("YES", "NO"))
#
ap_Inoculant.quantity$Inoculant.quantity <- droplevels(factor(ap_Inoculant.quantity$Inoculant.quantity))
# The number of Observations and StudyID
group_summary <- ap_Inoculant.quantity %>%
  group_by(Inoculant.quantity) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)
# Inoculant.quantity Observations Unique_StudyID
# <fct>                     <int>          <int>
# 1 NO                          861            154
# 2 YES                         236             57

overall_model_ap_Inoculant.quantity <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + inoculant.quantity, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Inoculant.quantity)
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127174.0291   254348.0582   254354.0582   254369.0999   254354.0798   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1531  0.3912    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1112) = 5587523.3541, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 612.9491, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# inoculant.quantityno     0.2898  0.0280  10.3612  <.0001   0.2350  0.3447  *** 
#   inoculant.quantityyes    0.0541  0.0289   1.8749  0.0608  -0.0025  0.1107    . 
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ap_Inoculant.quantity)
vcov_rotation <- vcov(overall_model_ap_Inoculant.quantity)
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
# Inoculant.quantityNO Inoculant.quantityYES 
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
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127438.8335   254877.6669   254889.6669   254919.7342   254889.7431   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1441  0.3796    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1109) = 2443191.9838, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 78.2654, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub      
# soil.origincropland      0.1740  0.0397  4.3785  <.0001   0.0961  0.2519  *** 
#   soil.originforest        0.1158  0.1167  0.9920  0.3212  -0.1130  0.3446      
# soil.origingrassland     0.2047  0.3801  0.5387  0.5901  -0.5402  0.9496      
# soil.originother         0.2635  0.0311  8.4649  <.0001   0.2025  0.3245  *** 
#   soil.originplantation    0.2625  0.1385  1.8951  0.0581  -0.0090  0.5339    . 
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
# soil.originother soil.originplantation  soil.origingrassland   soil.origincropland     soil.originforest 
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
#   1 not sandy                   110             25
# 2 sandy                       272             41

overall_model_ap_Sandy.or.not.sandy <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + sandy.or.not.sandy, random = ~ 1 | study.id, data = ap_Sandy.or.not.sandy, method = "REML")
# QM and p value
summary(overall_model_ap_Sandy.or.not.sandy)
# Multivariate Meta-Analysis Model (k = 382; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -14353.9266   28707.8533   28713.8533   28725.6738   28713.9171   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1572  0.3965     64     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 380) = 77742.8483, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 157.4886, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# sandy.or.not.sandynot sandy    0.3302  0.0513  6.4313  <.0001  0.2296  0.4308  *** 
#   sandy.or.not.sandysandy        0.1823  0.0510  3.5770  0.0003  0.0824  0.2822  *** 
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
# sandy.or.not.sandynot sandy     sandy.or.not.sandysandy 
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
# # A tibble: 2 × 3
# extraction.part Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 bulk                     643            120
# 2 rhizosphere              468             93

overall_model_ap_extraction.part <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + extraction.part, random = ~ 1 | study.id, data = ap_extraction.part, method = "REML")
# QM and p value
summary(overall_model_ap_extraction.part)
# Multivariate Meta-Analysis Model (k = 1111; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127379.7278   254759.4556   254765.4556   254780.4893   254765.4773   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1371  0.3703    204     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1109) = 2430301.5536, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 211.0810, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# extraction.partbulk           0.2022  0.0266   7.6046  <.0001  0.1501  0.2543  *** 
#   extraction.partrhizosphere    0.2690  0.0267  10.0852  <.0001  0.2167  0.3213  *** 
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
# extraction.partrhizosphere        extraction.partbulk 
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
# sterilization   Observations Unique_StudyID
# <fct>                  <int>          <int>
#   1 sterilization            397             61
# 2 unsterilization          717            151

overall_model_ap_Sterilization <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + sterilization, random = ~ 1 | study.id, data = ap_Sterilization, method = "REML")
# QM and p value
summary(overall_model_ap_Sterilization)

# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127437.4468   254874.8936   254880.8936   254895.9354   254880.9153   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1427  0.3777    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1112) = 5538878.9341, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 90.9357, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# sterilizationsterilization      0.2103  0.0273  7.6999  <.0001  0.1568  0.2639  *** 
#   sterilizationunsterilization    0.2376  0.0270  8.8066  <.0001  0.1847  0.2905  *** 
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
# sterilizationunsterilization   sterilizationsterilization 
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
# # A tibble: 2 × 3
# Stress    Observations Unique_StudyID
# <fct>            <int>          <int>
#   1 No Stress          938            201
# 2 Stress             167             29

overall_model_ap_Stress <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + stress, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Stress)
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -52730.2767  105460.5533  105466.5533  105481.5951  105466.5750   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1520  0.3899    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1112) = 3969090.2770, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 149505.5968, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# stressno stress    0.1904  0.0278   6.8567  <.0001  0.1360  0.2448  *** 
#   stressstress       0.5897  0.0278  21.2254  <.0001  0.5352  0.6441  *** 
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
# stressstress stressno stress 
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
#   1 native             133             27
# 2 no native          924            169

overall_model_ap_Native <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + native, random = ~ 1 | study.id, data = ap_Native, method = "REML")
# QM and p value
summary(overall_model_ap_Native)
# Multivariate Meta-Analysis Model (k = 1057; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -126953.7213   253907.4427   253913.4427   253928.3266   253913.4655   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1445  0.3801    190     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1055) = 5527417.9989, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 74.5762, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# nativenative       0.2182  0.0285  7.6478  <.0001  0.1623  0.2742  *** 
#   nativeno native    0.2334  0.0280  8.3240  <.0001  0.1785  0.2884  *** 
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
# nativeno native    nativenative 
# "a"             "b"











### 8.10 Plant.type
ap_Plant.type <- subset(ap, plant.type %in% c("crop", "herbaceous plants", "lianas", "shrub", "tree"))
#
ap_Plant.type$Plant.type <- droplevels(factor(ap_Plant.type$Plant.type))
# The number of Observations and StudyID
group_summary <-ap_Plant.type %>%
  group_by(Plant.type) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.ID)  
  )
print(group_summary)
# Plant.type        Observations Unique_StudyID
# <fct>                    <int>          <int>
#   1 crop                       365             67
# 2 herbaceous plants          553            103
# 3 lianas                      16              4
# 4 shrub                       76             17
# 5 tree                       104             18

overall_model_ap_Plant.type <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + plant.type, random = ~ 1 | study.id, data = ap, method = "REML")
# QM and p value
summary(overall_model_ap_Plant.type)
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127060.2081   254120.4162   254132.4162   254162.4835   254132.4925   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1613  0.4017    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1109) = 5560352.5591, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
#   QM(df = 5) = 833.3558, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub      
# plant.typecrop                 0.3905  0.0323  12.0780  <.0001   0.3271  0.4538  *** 
#   plant.typeherbaceous plants    0.1270  0.0320   3.9643  <.0001   0.0642  0.1898  *** 
#   plant.typelianas               0.2656  0.2024   1.3121  0.1895  -0.1312  0.6624      
# plant.typeshrub                0.3668  0.1019   3.6013  0.0003   0.1672  0.5665  *** 
#   plant.typetree                 0.0997  0.0930   1.0719  0.2838  -0.0826  0.2820      
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
# plant.typecrop             plant.typeshrub            plant.typelianas plant.typeherbaceous plants              plant.typetree 
# "a"                        "ab"                       "ab"                         "b"                        "b" 










### 8.11 Experimental.periods
ap_experimental.periods <- subset(ap, experimental.periods %in% c("Short-term", "Medium-term", "Long-term"))
#
ap_experimental.periods$experimental.periods <- droplevels(factor(ap_experimental.periods$experimental.periods))
# The number of Observations and StudyID
group_summary <-ap_experimental.periods %>%
  group_by(experimental.periods) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(study.id)  
  )
print(group_summary)
# A tibble: 3 × 3
# experimental.periods Observations Unique_StudyID
# <fct>                       <int>          <int>
#   1 Long-term                      55             14
# 2 Medium-term                   257             59
# 3 Short-term                    622            113

overall_model_Experimental.periods <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + experimental.periods, random = ~ 1 | study.id, data = ap_experimental.periods, method = "REML")
# QM and p value
summary(overall_model_Experimental.periods)
# Multivariate Meta-Analysis Model (k = 934; method: REML)
# 
# logLik     Deviance          AIC          BIC         AICc   
# -98728.2347  197456.4693  197464.4693  197483.8144  197464.5125   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1677  0.4095    177     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 931) = 821663.3072, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
#   QM(df = 3) = 4457.0236, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# experimental.periodsLong-term      0.7298  0.0345  21.1514  <.0001  0.6622  0.7975  *** 
#   experimental.periodsMedium-term    0.3054  0.0340   8.9840  <.0001  0.2388  0.3721  *** 
#   experimental.periodsShort-term     0.1241  0.0323   3.8467  0.0001  0.0609  0.1873  *** 
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
# experimental.periodslong-term experimental.periodsmedium-term  experimental.periodsshort-term 
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
# Multivariate Meta-Analysis Model (k = 1114; method: REML)
# 
# logLik      Deviance           AIC           BIC          AICc   
# -127392.7876   254785.5753   254791.5753   254806.6170   254791.5969   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2    0.1422  0.3771    205     no  study.id 
# 
# Test for Residual Heterogeneity:
#   QE(df = 1112) = 5287744.7939, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
#   QM(df = 2) = 186.7924, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# fertilizerfertilizer       0.2323  0.0269  8.6416  <.0001  0.1796  0.2850  *** 
#   fertilizerno fertilizer    0.2292  0.0269  8.5255  <.0001  0.1765  0.2818  *** 
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
# fertilizerfertilizer fertilizerno fertilizer 
# "a"                     "b" 
library(patchwork)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+ plot_layout(ncol = 3)
#图片导出为23×20






library(lme4)
ap$Wr <- 1 / ap$Vi
###background character
##time
sum(!is.na(ap$experimental.duration))#934

ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + inoculant.type + (1|study.id), weights=Wr, data=ap)
m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + (1|study.id), weights=Wr, data=ap)
anova(m1, m2) 
# m1: RR ~ dur_z + inoculant.type + (1 | study.id)
# m2: RR ~ dur_z + I(dur_z^2) + inoculant.type + (1 | study.id)
# npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
# m1    7 2585.3 2619.2 -1285.6    2571.3                         
# m2    8 2571.3 2610.0 -1277.7    2555.3 15.968  1  6.442e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z          5027.1  5027.1     1 482.18 20.5668 7.275e-06 ***
# I(dur_z^2)     3937.2  3937.2     1 622.70 16.1075 6.713e-05 ***
# inoculant.type 3266.0  1088.7     3 336.93  4.4539  0.004378 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b1 <- fixef(m2)["dur_z"]
b2 <- fixef(m2)["I(dur_z^2)"]

z_star <- -b1/(2*b2)

mu  <- attr(scale(ap$experimental.duration), "scaled:center")
sdv <- attr(scale(ap$experimental.duration), "scaled:scale")

t_star <- mu + z_star*sdv
t_star
# dur_z 
# 1513.208 
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
sum(!is.na(ap$experimental.duration)) #934
# 图片导出为8*8









#pH
ap$pH_z <- as.numeric(scale(ap$soil.ph))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + (1|study.id),
           weights = Wr, data = ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + I(pH_z^2) + (1|study.id),
           weights = Wr, data = ap)

anova(m1, m2)   # 看二次项是否显著改善拟合
# refitting model(s) with ML (instead of REML)
# Data: ap
# Models:
#   m1: RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + (1 | study.id)
# m2: RR ~ dur_z + I(dur_z^2) + inoculant.type + pH_z + I(pH_z^2) + (1 | study.id)
# npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)    
# m1    9 1221.5 1262.8 -601.74   1203.48                         
# m2   10 1006.8 1052.7 -493.41    986.82 216.66  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# dur_z              0.3     0.3     1 264.57   0.0062    0.9372    
# I(dur_z^2)        35.3    35.3     1 551.13   0.8084    0.3690    
# inoculant.type  1269.5   423.2     3 683.15   9.6938 2.852e-06 ***
# pH_z           22761.6 22761.6     1 686.32 521.4258 < 2.2e-16 ***
# I(pH_z^2)      11352.2 11352.2     1 692.86 260.0572 < 2.2e-16 ***
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
            linewidth = 1.2) +
  geom_vline(xintercept = pH_star, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  theme_bw() +
  labs(x = "Soil pH", y = "Effect size (RR of Soil-available P)")
p16
sum(!is.na(ap$pH))#869







##TC
library(ggeffects)
library(ggplot2)
ap$TC_z <- as.numeric(scale(ap$soil.tc))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + TC_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + TC_z + I(TC_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)   # 若不显著，停止
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# dur_z            3.318   3.318     1    39  0.4898 0.4881658    
# I(dur_z^2)      24.094  24.094     1    39  3.5569 0.0667608 .  
# inoculant.type  36.277  18.138     2    39  2.6777 0.0813427 .  
# TC_z           116.207 116.207     1    39 17.1550 0.0001788 ***
# I(TC_z^2)      152.471 152.471     1    39 22.5085 2.797e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 3. 计算拐点
b1 <- fixef(m2)["TC_z"]
b2 <- fixef(m2)["I(TC_z^2)"]

z_star  <- -b1 / (2 * b2)

mu  <- mean(ap$soil.tc, na.rm = TRUE)
sdv <- sd(ap$soil.tc, na.rm = TRUE)

tc_star <- z_star * sdv + mu
tc_star

## 4. 生成平滑预测曲线
ap_used <- model.frame(m2)
zmin <- min(ap_used$TC_z, na.rm = TRUE)
zmax <- max(ap_used$TC_z, na.rm = TRUE)

pred_z <- ggpredict(
  m2,
  terms = sprintf("TC_z [%.3f:%.3f by=0.01]", zmin, zmax)
)

## 5. 把 z 值换回原始 soil.tc
pred_z$soil.tc <- pred_z$x * sdv + mu

## 6. 提取模型实际使用的数据行来画散点
idx <- as.integer(rownames(ap_used))
ap_used$soil.tc <- ap$soil.tc[idx]

## 7. 作图
p17 <- ggplot() +
  geom_point(
    data = ap_used,
    aes(x = soil.tc, y = RR),
    color = "gray",
    size = 10,
    shape = 21
  ) +
  geom_ribbon(
    data = pred_z,
    aes(x = soil.tc, ymin = conf.low, ymax = conf.high),
    alpha = 0.20
  ) +
  geom_line(
    data = pred_z,
    aes(x = soil.tc, y = predicted),
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = tc_star,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  theme_bw() +
  labs(
    x = "Soil TC",
    y = "Effect size (RR of Soil-available P)"
  )

p17










##SOC

ap$soc_z <- as.numeric(scale(ap$soil.soc))
ap$dur_z <- as.numeric(scale(ap$experimental.duration))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + soc_z + I(soc_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)   # 若不显著，停止
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            9.78   9.778     1 160.76  0.3472   0.55653    
# I(dur_z^2)     127.81 127.808     1 211.04  4.5380   0.03431 *  
# inoculant.type 666.36 222.120     3 223.61  7.8866 5.036e-05 ***
# soc_z          156.83 156.834     1 210.93  5.5686   0.01920 *  
# I(soc_z^2)     183.11 183.109     1 135.76  6.5015   0.01189 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b1 <- fixef(m2)["soc_z"]
b2 <- fixef(m2)["I(soc_z^2)"]
z_star <- -b1/(2*b2)

mu  <- mean(ap$soil.soc, na.rm=TRUE)
sdv <- sd(ap$soil.soc, na.rm=TRUE)
soc_star <- mu + z_star*sdv

c(z_star=z_star, soc_star=soc_star)
library(ggeffects)
library(ggplot2)

pred <- ggpredict(m2, terms = "soc_z [all]")

mu  <- mean(ap$soil.soc, na.rm=TRUE)
sdv <- sd(ap$soil.soc, na.rm=TRUE)
pred$SOC <- pred$x * sdv + mu

# 拐点
b1 <- fixef(m2)["soc_z"]
b2 <- fixef(m2)["I(soc_z^2)"]
z_star <- -b1/(2*b2)
soc_star <- mu + z_star*sdv
soc_star
# 用模型实际使用的行画点（避免 NA）
ap_used <- model.frame(m2)
idx <- as.integer(rownames(ap_used))
ap_used$soil.soc <- ap$soil.soc[idx]

p18 <- ggplot() +
  geom_point(data=ap_used, aes(x=soil.soc, y=RR),
             color="gray", size=10, shape = 21) +
  geom_ribbon(data=pred, aes(x=SOC, ymin=conf.low, ymax=conf.high), alpha=0.2) +
  geom_line(data=pred, aes(x=SOC, y=predicted), linewidth=1.2) +
  geom_vline(xintercept=soc_star, linetype="dashed", linewidth=0.8) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.6) +
  theme_bw() +
  labs(x="Soil SOC", y="Effect size (RR / lnRR of available P)")

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
#                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z           73.632  73.632     1 54.587  4.9730   0.02988 *  
# I(dur_z^2)     312.771 312.771     1 47.338 21.1239 3.212e-05 ***
# inoculant.type 115.482  38.494     3 97.464  2.5998   0.05651 .  
# tn_z             3.144   3.144     1 40.220  0.2124   0.64741    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            0.464   0.464     1 38.293  0.0497 0.8247449    
# I(dur_z^2)      19.310  19.310     1 43.032  2.0685 0.1576062    
# inoculant.type 183.346  61.115     3 57.584  6.5465 0.0006933 ***
# tp_z             3.487   3.487     1 20.074  0.3735 0.5479479    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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







##TK
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)

# 标准化
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$tk_z  <- as.numeric(scale(ap$tk))

# 模型比较
m1 <- lmer(
  RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + (1 | study.id),
  weights = Wr, data = ap
)

m2 <- lmer(
  RR ~ dur_z + I(dur_z^2) + inoculant.type + tk_z + I(tk_z^2) + (1 | study.id),
  weights = Wr, data = ap
)

anova(m1, m2)
anova(m2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)  
# dur_z           40.564  40.564     1  5.7617  2.6439 0.1571  
# I(dur_z^2)      17.405  17.405     1  5.4158  1.1345 0.3320  
# inoculant.type 181.148  60.383     3 11.7214  3.9357 0.0370 *
# tk_z            20.653  20.653     1  9.9263  1.3461 0.2731  
# I(tk_z^2)       16.905  16.905     1  8.0847  1.1019 0.3242  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# -----------------------------
# 计算 tk 的二次曲线拐点
# -----------------------------
b1 <- fixef(m2)["tk_z"]
b2 <- fixef(m2)["I(tk_z^2)"]

z_star <- -b1 / (2 * b2)

mu  <- mean(ap$tk, na.rm = TRUE)
sdv <- sd(ap$tk, na.rm = TRUE)

tk_star <- mu + z_star * sdv

c(z_star = z_star, tk_star = tk_star)

# -----------------------------
# 生成预测值
# -----------------------------
pred <- ggpredict(m2, terms = "tk_z [all]")

# 把标准化坐标换回原始 tk 尺度
pred$tk <- pred$x * sdv + mu

# -----------------------------
# 取模型实际使用的数据行（避免 NA 干扰）
# -----------------------------
ap_used <- model.frame(m2)
idx <- as.integer(rownames(ap_used))
ap_used$tk <- ap$tk[idx]

# -----------------------------
# 作图
# -----------------------------
p21 <- ggplot() +
  geom_point(data=ap_used, aes(x=tk, y=RR),
             color="gray", size=10, shape=21) +
  geom_ribbon(data=pred, aes(x=tk, ymin=conf.low, ymax=conf.high), alpha=0.2) +
  geom_line(data=pred, aes(x=tk, y=predicted), linewidth=1.2) +
  geom_vline(xintercept=tk_star, linetype="dashed", linewidth=0.8) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.6) +
  theme_bw() +
  labs(x="Soil TK", y="Effect size (RR / lnRR of available P)")

p21











##AN
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$an_z  <- as.numeric(scale(ap$an))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + an_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + an_z + I(an_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# dur_z           55.62  55.618     1  3.882  1.7018 0.264021   
# I(dur_z^2)       3.44   3.444     1  7.783  0.1054 0.754021   
# inoculant.type 629.13 209.711     3 42.706  6.4167 0.001098 **
# an_z             8.83   8.828     1  4.306  0.2701 0.628843   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p22 <- ggplot(ap, aes(y=RR, x=an)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil AN" , y="RR of Soil-available P")
p22





##AP
ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ap_z  <- as.numeric(scale(ap$ap))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ap_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~dur_z + I(dur_z^2) + inoculant.type + ap_z + I(ap_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)


anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# dur_z            2.99    2.99     1 162.18  0.1248    0.7244    
# I(dur_z^2)      15.20   15.20     1 312.97  0.6334    0.4267    
# inoculant.type 957.09  319.03     3 361.41 13.2987 2.974e-08 ***
# ap_z             0.00    0.00     1  58.25  0.0001    0.9931    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p23 <- ggplot(ap, aes(y=RR, x=ap)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil AP" , y="RR of Soil-available P")
p23






##AK

ap$dur_z <- as.numeric(scale(ap$experimental.duration))
ap$ak_z  <- as.numeric(scale(ap$ak))

m1 <- lmer(RR ~ dur_z + I(dur_z^2) + inoculant.type + ak_z + (1|study.id),
           weights=Wr, data=ap)

m2 <- lmer(RR ~ dur_z + inoculant.type + ak_z + I(ak_z^2) + (1|study.id),
           weights=Wr, data=ap)

anova(m1, m2)
anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# dur_z          35.526  35.526     1  53.185  4.4688 0.03923 *
# I(dur_z^2)     14.270  14.270     1  47.862  1.7951 0.18663  
# inoculant.type 22.468   7.489     3  84.481  0.9421 0.42410  
# ak_z           29.673  29.673     1 147.994  3.7325 0.05527 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p24 <- ggplot(ap, aes(y=RR, x=ak)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Soil AK" , y="RR of Soil-available P")
p24
#图片输出为11×12
p16+p17+p18+p19+p20+p21+p22+p23+p24





# Unit: cells/seed
Bacteria_cell_seed <-read.csv(file.choose())
Bacteria_cell_seed$Wr <- 1 / Bacteria_cell_seed$Vi
mdBacteria_cell_seed<-lmer(RR~ scale(final.inoculation.amount)+scale(experimental.duration) + (1|study.id), weights=Wr, data=Bacteria_cell_seed)
anova(mdBacteria_cell_seed) 
r.squaredGLMM(mdFinalinoculationquantification)
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)
r.squaredGLMM(mdLog10)








# Unit: Bacteria_cfu_ha_soil
Bacteria_cfu_ha_soil <-read.csv(file.choose())
Bacteria_cfu_ha_soil$Wr <- 1 / Bacteria_cfu_ha_soil$Vi
##Finalinoculationquantification
mdBacteria_cfu_ha_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_ha_soil)
anova(mdBacteria_cfu_ha_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(experimental.duration)    1.26635 1.26635     1 4.076  4.9161 0.08963 .
# scale(final.inoculation.amount) 0.00244 0.00244     1 1.014  0.0095 0.93806  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

p25 <- ggplot(Bacteria_cfu_ha_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="Cfu/ha soil" , y="RR of Soil-available P")
p25
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.003183279 0.3380941
##Log10
mdlBacteria_cfu_ha_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_ha_soil)
anova(mdlBacteria_cfu_ha_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(experimental.duration) 1.26635 1.26635     1 4.076  4.9161 0.08963 .
# scale(log10)                 0.00244 0.00244     1 1.014  0.0095 0.93806  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p26 <- ggplot(Bacteria_cfu_ha_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  labs(x="log10(Cfu/ha soil)" , y="RR of Soil-available P")
p26







# Unit: Bacteria_cfu_kg_soil
Bacteria_cfu_kg_soil <-read.csv(file.choose())
Bacteria_cfu_kg_soil$Wr <- 1 / Bacteria_cfu_kg_soil$Vi
##Finalinoculationquantification
mdBacteria_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_kg_soil)
anova(mdBacteria_cfu_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(experimental.duration)    7.6042  7.6042     1 25.104  1.7121 0.2026
# scale(final.inoculation.amount) 2.4331  2.4331     1 10.061  0.5478 0.4761
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
#                               Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(experimental.duration) 1.26635 1.26635     1 4.076  4.9161 0.08963 .
# scale(log10)                 0.00244 0.00244     1 1.014  0.0095 0.93806  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p28 <- ggplot(Bacteria_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="log10(Cfu/kg soil)" , y="RR of Soil-available P")
p28





# Unit: Bacteria_cfu_plant
Bacteria_cfu_plant <-read.csv(file.choose())
Bacteria_cfu_plant$Wr <- 1 / Bacteria_cfu_plant$Vi
##Finalinoculationquantification
mdBacteria_cfu_plant<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Bacteria_cfu_plant)
anova(mdBacteria_cfu_plant)
r.squaredGLMM(mdFinalinoculationquantification)
mdlBacteria_cfu_plant<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Bacteria_cfu_plant)
anova(mdlBacteria_cfu_kg_soil)







# Unit: Fungi_cfu_kg_soil
Fungi_cfu_kg_soil <-read.csv(file.choose())
Fungi_cfu_kg_soil$Wr <- 1 / Fungi_cfu_kg_soil$Vi
##Finalinoculationquantification
mdFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(final.inoculation.amount)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdFungi_cfu_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(experimental.duration)    0.70721 0.70721     1 3.2807  0.0214 0.8922
# scale(final.inoculation.amount) 0.54521 0.54521     1 3.1780  0.0165 0.9055
p29 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  theme(panel.grid=element_blank())+ 
  labs(x="cfu/kg soil" , y="RR of Soil-available P")
p29
r.squaredGLMM(mdFinalinoculationquantification)
mdlFungi_cfu_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_cfu_kg_soil)
anova(mdlFungi_cfu_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(experimental.duration)  31.768  31.768     1     8  1.0423 0.33717  
# scale(log10)                 114.036 114.036     1     8  3.7413 0.08913 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
p30 <- ggplot(Fungi_cfu_kg_soil, aes(y=RR, x=log10)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
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
mdlFungi_pieces_root<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=Fungi_pieces_root)
anova(mdlFungi_pieces_root)








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
#                                 Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration)    24.046  24.046     1     2  4.7047 0.1623
# scale(final.inoculation.amount)  4.125   4.125     1     2  0.8071 0.4638
p31 <- ggplot(AMF_propagules_kg_soil, aes(y=RR, x=final.inoculation.amount)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="solid") + 
  theme_bw()+
  theme(text = element_text(                                                                                                                                                                              size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="solid")+
  theme(panel.grid=element_blank())+ 
  labs(x="propagules/kg soil" , y="RR of Soil-available P")
p31

r.squaredGLMM(mdFinalinoculationquantification)
mdlAMF_propagules_kg_soil<-lmer(RR~ scale(experimental.duration)+scale(log10)+ (1|study.id), weights=Wr, data=AMF_propagules_kg_soil)
anova(mdlAMF_propagules_kg_soil)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(experimental.duration) 24.046  24.046     1     2  4.7047 0.1623
# scale(log10)                  4.125   4.125     1     2  0.8071 0.4638
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
# scale(experimental.duration)    5.0326  5.0326     1  5.524  0.1827 0.6852
# scale(final.inoculation.amount) 2.5244  2.5244     1 42.687  0.0916 0.7636
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
#                               Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(experimental.duration)  1.7906  1.7906     1  5.3322  0.0652 0.8080
# scale(log10)                 19.8313 19.8313     1 16.9613  0.7226 0.4071
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

# 如果你更喜欢“每个变量一页 PDF”，用这个：
# pdf("RR_vs_traits_HolmAdjusted_multipage.pdf", width = 8, height = 8)
# for (p in plots) print(p)
# dev.off()


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
# soil.ph               60.371990   0.00990099   68.04124653         0.00990099
# experimental.duration 58.306042   0.00990099   49.28984029         0.00990099
# soil.soc              40.909463   0.00990099   20.49939988         0.00990099
# ap                    39.309446   0.00990099   16.52709284         0.00990099
# tp                    23.051993   0.00990099    5.67772804         0.32673267
# tn                    21.211299   0.00990099   11.06128477         0.00990099
# ak                    18.651140   0.00990099    5.58573883         0.44554455
# an                    12.286023   0.00990099    1.29646016         0.98019802
# tk                    11.835977   0.00990099    1.77676656         0.80198020
# soil.tc                8.849218   0.03960396    1.32337450         0.54455446
# organic.p              0.000000   1.00000000    0.02426327         0.80198020
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
  limits = c(0, 65),
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
#       30      311       52 
## 8*8
