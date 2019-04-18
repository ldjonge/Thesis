#!/usr/bin/Rscript
setwd("C:/Users/Lorenzo/Documents/Thesis/scripts")

table = read.csv("outTable", sep="\t", header=TRUE)
par(mfrow=c(1,2))
plot(1:nrow(table), table$A, type="l", col="blue", ylim=c(0,1), xlab="Generation", ylab="Frequency")
lines(1:nrow(table), table$I, col="green")
lines(1:nrow(table), table$O, col = "red")
legend("topright", legend=c("A", "I", "O"), col=c("blue", "green", "red"), lty=1)

fullPop = table$Males + table$Females
plot(fullPop[2:nrow(table)], xlab="Generation", ylab="N", ylim=c(0,300), type="l", lty=8)
lines(1:nrow(table), table$Males, type = "l", col="blue")
lines(1:nrow(table), table$Females, col="red")
legend("topright", legend=c("Male", "Female", "Total"), col=c("blue", "red", "black"), lty=c(1,1,8))
