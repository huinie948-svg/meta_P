library(ggplot2)
library(maps)
library(dplyr)

world.dat <- map_data("world")

world.map <- ggplot() +
  geom_polygon(data = world.dat,
               aes(x = long, y = lat, group = group),
               fill = "#dedede", color = "white", linewidth = 0.2) +
  coord_quickmap() +
  theme_void()

df <- read.csv(file.choose(), stringsAsFactors = FALSE)

# 清理经纬度，避免负号读错
df$Longitude <- as.numeric(gsub("[−–]", "-", df$Longitude))
df$Latitude  <- as.numeric(gsub("[−–]", "-", df$Latitude))

# 如果你只有一种类型，直接指定一个常量列（用于图例/颜色）
df$EWEs <- "Cold"

p <- world.map +
  geom_jitter(
    data = df,
    aes(x = Longitude, y = Latitude, fill = EWEs),
    width = 0.15, height = 0.15,     # 抖动幅度（单位：度）
    shape = 1,
    size = 4,
    color = "#448DCD",
    stroke = 0.25,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("Cold" = "#448DCD"), name = "EWEs") +
  theme(
    legend.position = c(0.10, 0.30),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )

p

