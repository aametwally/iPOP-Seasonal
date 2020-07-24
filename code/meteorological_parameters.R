library(ggplot2)
library(lubridate)
library(reshape2)


weath_df = read.csv("meteorological_parameters.csv")
weath_df$day = as.Date(weath_df$day, format = '%m/%d/%y')
weath_df$dayofyear = yday(weath_df$day)
weath_df_selected = weath_df[,c("dayofyear", "avgTemp", "avgHumid", "avgSolarRad", 
                                "avgPress", "avgWind", "avgPrecipMidNight")]
weath_df_selected$avgPress = weath_df_selected$avgPress/1000

df = melt(weath_df_selected, id.vars = "dayofyear")
variable.labs <- c("Temperature (°C)", "Humidity (%)", "Solar Rad.(W/m^2)", "Pressure (kPa)", "Wind Speed (m/s)", "Precipitation (mm)")
names(variable.labs) <- c("avgTemp", "avgHumid", "avgSolarRad", "avgPress", "avgWind", "avgPrecipMidNight")


pdf("meteorological_parameters.pdf", h = 12, w = 8)
ggplot(df, aes(dayofyear, value)) + 
  geom_point(stat='identity', size = 1, show.legend = TRUE, alpha = 0.6, color = "black") +
  facet_grid(rows = vars(variable), scales = "free_y", labeller = labeller(variable = variable.labs)) +
  theme_bw() + 
  geom_smooth() +
  labs(y = "", x = "Day of the Year") +
  ggtitle("Meteorological Parameters") + 
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size=15, face="bold"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"))
dev.off()