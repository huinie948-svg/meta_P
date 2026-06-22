library(readr)
library(dplyr)
library(janitor)

# 1) 读取两张表（把路径改成你的实际文件名/路径）
t1 <- read_csv("ap.csv", show_col_types = FALSE) %>%
  clean_names()

t2 <- read_csv("phosphatase.csv", show_col_types = FALSE) %>%
  clean_names()

# 2) 计算 study_id 交集（注意：clean_names() 后 study id 会变成 study_id）
id_inter <- intersect(
  unique(na.omit(t1$study_id)),
  unique(na.omit(t2$study_id))
)

# 3) 过滤到交集
t1_i <- t1 %>%
  filter(study_id %in% id_inter) %>%
  mutate(source_table = "table1")

t2_i <- t2 %>%
  filter(study_id %in% id_inter) %>%
  mutate(source_table = "table2")

merged <- bind_rows(t1_i, t2_i)
front <- c("study_id", "source_table")
merged <- merged %>% select(any_of(front), everything())

# 6) 输出
write_csv(merged, "merged_intersection_by_study_id.csv")


################################################  piecewiseSEM
sem_data <- read.csv("data_collection.csv", fileEncoding = "latin1")
library(piecewiseSEM)
library(nlme)
library(dplyr)

# 读入
sem_data <- read.csv("data_collection.csv", fileEncoding = "latin1")
names(sem_data)
# 选列（你的写法）
sem_data <- sem_data[, c("ï..ID", "soil.ph", "lnphosphatase", "lnap",
                         "RRbiomass", "RR_pplant", "RR_proot","RRshoot_biomass","RRroot_biomass"
                         ,"RRleaf_biomass", "RRyield"
                         )]

# 重命名 ID，避免乱码列名后面一直麻烦
names(sem_data)[1] <- "ID"

# ✅ 关键：ID 保持为因子；其它列转 numeric
sem_data <- sem_data %>%
  mutate(
    ID = as.factor(ID),
    across(-ID, ~ as.numeric(as.character(.)))
  )

sem_data$ph_z  <- scale(sem_data$soil.ph)
sem_data$ph_z2 <- sem_data$ph_z^2


# 子模型 1：pH -> AP
m_ap <- lme(
  lnap ~ lnphosphatase+ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)

# 子模型 2：pH -> phosphatase
m_phos <- lme(
  lnphosphatase ~ ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)

# 子模型 3：AP + phosphatase -> biomass
m_bio <- lme(
  RRbiomass ~ lnap + lnphosphatase,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)
sem_mod <- psem(m_ap, m_phos, m_bio)

summary(sem_mod)                        # Fisher's C、路径显著性等
coefs(sem_mod, standardize = "scale")   # 标准化路径系数（推荐看这个）
rsquared(sem_mod)                       # 每个子模型的 marginal/conditional R2


# 子模型 1：pH -> AP
m1_ap <- lme(
  lnap ~ lnphosphatase+ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)

# 子模型 2：pH -> phosphatase
m1_phos <- lme(
  lnphosphatase ~ ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)

# 子模型 3：AP + phosphatase -> biomass
m1_bio <- lme(
  RRroot_biomass ~ lnap + lnphosphatase,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit
)
sem1_mod <- psem(m1_ap, m1_phos, m1_bio)

summary(sem1_mod)                        # Fisher's C、路径显著性等
coefs(sem1_mod, standardize = "scale")   # 标准化路径系数（推荐看这个）
rsquared(sem1_mod)                       # 每个子模型的 marginal/conditional R2


# 子模型 1
m3_ap <- lme(
  lnap ~ lnphosphatase+ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 2
m3_phos <- lme(
  lnphosphatase ~ ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 3
m3_bio <- lme(
  RRshoot_biomass ~ lnap + lnphosphatase,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)
sem3_mod <- psem(m3_ap, m3_phos, m3_bio)
summary(sem3_mod)                        # Fisher's C、路径显著性等
coefs(sem3_mod, standardize = "scale")   # 标准化路径系数（推荐看这个）
rsquared(sem3_mod)   



# 子模型 1
m4_ap <- lme(
  lnap ~ lnphosphatase+ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 2
m4_phos <- lme(
  lnphosphatase ~ ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 3
m4_bio <- lme(
  RRyield ~ lnap + lnphosphatase+ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)
sem4_mod <- psem(m4_ap, m4_phos, m4_bio)
summary(sem4_mod)                        # Fisher's C、路径显著性等
coefs(sem4_mod, standardize = "scale")   # 标准化路径系数（推荐看这个）
rsquared(sem4_mod)  
names(sem_data)



# 子模型 1
m5_ap <- lme(
  lnap ~ lnphosphatase+ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 2
m5_phos <- lme(
  lnphosphatase ~ ph_z + ph_z2,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)

# 子模型 3
m5_bio <- lme(
  RRleaf_biomass ~ lnap + lnphosphatase,
  random = ~ 1 | ID,
  data = sem_data,
  method = "REML",
  na.action = na.omit)
sem5_mod <- psem(m5_ap, m5_phos, m5_bio)
summary(sem5_mod)                        # Fisher's C、路径显著性等
coefs(sem5_mod, standardize = "scale")   # 标准化路径系数（推荐看这个）
rsquared(sem5_mod)  
names(sem_data)









