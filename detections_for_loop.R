## Counting seconds
setwd("C:/Users/Karin/Desktop")
data <- read.csv("sgl1_20170604_063001.detectoR.csv", h = T)

detections <- NULL

for (i in 1:length(data[, 1])) {
  
  if (i == length(data[, 1])) {
    
    detections[i] <- paste("END")
  
  } else if (abs(data$Begin.Time[i] - data$Begin.Time[i + 1]) < 1) {
    
    detections[i] <- paste("Good")
  
  } else if (abs(data$Begin.Time[i] - data$Begin.Time[i + 1]) >= 1) {
    
    detections[i] <- paste("Bad")
 }
}

data$detections <- detections
