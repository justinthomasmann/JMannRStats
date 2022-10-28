## warbleR test 1

require(warbleR)
?imp.raven()
setwd("C:/Users/Karin/Desktop/Recordings")
checkwavs(path = NULL)


manualoc()

?specan()

#### example

setwd(tempdir())
data("selection.files")
write.table(selection.files[[1]],file = "100889-Garrulax monileger.selections.txt",
            row.names = FALSE, sep= "\t")

write.table(selection.files[[2]],file = "1023-Arremonops rufivirgatus.selections.txt",
            row.names = FALSE, sep= "\t")

#providing the name of the column with the sound file names
rav.dat<-imp.raven(sound.file.col = "End.File", all.data = FALSE)

View(rav.dat)

#getting all the data
rav.dat<-imp.raven(all.data = TRUE)
View(rav.dat)



in

?specan()

data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", "selec.table"))
writeWave(Phae.long1,"Phae.long1.wav")
writeWave(Phae.long2,"Phae.long2.wav")
writeWave(Phae.long3,"Phae.long3.wav")
writeWave(Phae.long4,"Phae.long4.wav")

a <- specan(X = selec.table, bp = c(0, 22))
a
# using a diferent threshold
a <- specan(X = selec.table, bp = c(0, 22), threshold = 20)
View(a)
