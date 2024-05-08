library(tidyverse)
library(MatchIt)
library(optmatch)

setwd("D:/Dissertation")

df <- read_csv('293B_ControlRandom.csv')


df_case <- df %>% 
  filter(Group == 1)

df_control <- df %>% 
  filter(Group == 0) %>% 
  anti_join(df_case)

final_df <- bind_rows(df_case, df_control)
final_df <- final_df[complete.cases(final_df$HS_GPA),]
final_df <- final_df[complete.cases(final_df$M_SAT),]

matches_out <- matchit(Group ~ Male + HS_GPA + M_SAT, 
                       data = final_df,
                       method = "optimal",
                       ratio = 2
)
matches_out
summary(matches_out)
plot(matches_out)

df_analysis <- match.data(matches_out)

write.csv(df_analysis,"D:\\Dissertation\\RandomControl.csv", row.names = TRUE)

