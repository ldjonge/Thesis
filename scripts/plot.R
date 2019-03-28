#!/usr/bin/Rscript
setwd("C:/Users/Lorenzo/Documents/Thesis/scripts")

table = read.csv("outTable", sep="\t", header=TRUE)

plot(1:nrow(table), table$A, type="l", col="blue", ylim=c(0,1), xlab="Generation", ylab="Frequency")
lines(1:nrow(table), table$I, col="green")
lines(1:nrow(table), table$O, col = "red")
legend("topright", legend=c("A", "I", "O"), col=c("blue", "green", "red"), lty=1)

plot(1:nrow(table), table$Males, type = "l", col="blue")
lines(1:nrow(table), table$Females, col="red")
fullPop = table$Males + table$Females
plot(fullPop[2:nrow(table)])
