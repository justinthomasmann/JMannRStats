setwd("C:/Users/Justin/Desktop/ARUData/truth_tables")

t.table <- read.csv("sgl6_20170530_073002.TruthTable.selections.csv", h=T)

singFreq.df <- data.frame(t.table$Selection, t.table$beginTime, difference)

difference <- diff(t.table$beginTime)
difference


singFreq.df$difference <- difference

singFreq.df[-1, ] - singFreq.df[-nrow(singFreq.df), ]
