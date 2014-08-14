###############################################################################################################
#this code produces a dynamic animation of the reaction-diffusion system simulated in the accompanying 
#Python code. to run this code, you need to generate the "stateX_n=100.csv" data using the Python code
#found in "Turing_patterns_OOP_implemantation.ipynb"
###############################################################################################################

library(animation)

library(data.table)
# Start the clock!
ptm <- proc.time()
setwd("~/Documents/")
mydata = fread("stateX_n=200.csv", header=FALSE, sep = ',')

# Stop the clock
proc.time() - ptm

mydata1 <- mydata



system.time(mydata = read.csv("stateY_n=150.csv", header=FALSE))# read csv file 


library(animation)

setwd("~/Desktop/")
mydata = read.csv("stateY_n=150.csv", header=FALSE)

my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(min(mydata),min(mydata)/2,length=100),  # for blue
               seq(min(mydata)/2,0,length=100),              # for yellow
               seq(0,max(mydata),length=100))              # for red


ani.record(reset = TRUE) # clear history before recording
# 

#oopt = ani.options(interval = 0.15, nmax = 200)

oopt = ani.options(interval = 0.15, nmax = length(mydata))
## use a loop to create images one by one
for (i in seq(1,ani.options("nmax"), 5)) {
  
  a <- matrix(mydata[[i]],nrow=250,byrow = TRUE)
  heatmap( a , Rowv = NA , Colv = NA, labRow = NA , labCol = NA , scale = "none", col=my_palette)
  ani.record() # record the current frame
  ani.pause() ## pause for a while (’interval’)
  print(i)
}

saveHTML(ani.replay() , outdir = getwd() , img.name = "record_plot")

saveGIF(ani.replay() , outdir = getwd() , img.name = "record_plot")
