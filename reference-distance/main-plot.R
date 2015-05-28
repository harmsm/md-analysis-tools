files <- c("../SR2-Q41E-M75L_N0600_run000-processed.txt",
           "../SR2-Q41E-M75L_N0600_run001-processed.txt",
           "../SR2-Q41E-M75L_N0600_run002-processed.txt",
           "../SR2_N0600_run000-processed.txt",
           "../SR2_N0600_run001-processed.txt",
           "../SR2_N0600_run002-processed.txt")
#          "../SR2-Q41E-M75L_N1500_run000-processed.txt",
#          "../SR2-Q41E-M75L_N1500_run001-processed.txt",
#          "../SR2-Q41E-M75L_N1500_run002-processed.txt",
#          "../SR2_N1500_run000-processed.txt",
#          "../SR2_N1500_run001-processed.txt",
#          "../SR2_N1500_run002-processed.txt")

x <- seq(0,1,length=length(files))
colors <- c("blue","slateblue","lightblue","red","orange","peachpuff")

# Only plot histograms for frames above this value (for removing burn in)
frame_cutoff <- 0

# Pool histograms for repeated trajectories (assumes you read them in in order
# and that there are three reps for each trajectory)
pool_histograms <- TRUE


d <- subset(read.table(files[1]),frame>frame_cutoff)
columns <- colnames(d)
columns <- columns[2:length(columns)]

num_snapshots <- length(d[,1])

# Go through every distance interaction
for (i in 1:length(columns)){

    print(columns[i])

    # Figure out the relevant breaks for the histogram based on the minium and
    # maximum distances observed over all trajectories
    out_matrix <- matrix(NA,length(files),num_snapshots)
    for (j in 1:length(files)){
        
        d <- subset(read.table(files[j]),frame>frame_cutoff)
        out_matrix[j,] <- d[,i+1]

    }

    lim <- c(min(out_matrix),max(out_matrix))
    breaks <- seq(lim[1]-5,lim[2]+5,by=0.5)

    # Start a new plot
    pdf(paste(columns[i],".pdf",sep=""))

    plot(NA,xlim=lim,ylim=c(0,0.5),axes=FALSE,ann=FALSE)
    axis(1,lwd=1.5); mtext("distance (A)",1,line=2.0,cex=1.2)
    axis(2,lwd=1.5); mtext("frequency",2,line=2.0,cex=1.2)
    mtext(columns[i],3,line=1.5,cex=1.5)

    # Go through every trajectory and plot a histogram for the given distance
    for (j in 1:length(files)){
        
        a <- hist(out_matrix[j,],plot=FALSE,breaks=breaks)

        # If we're pooling replicate trajectories...
        if (pool_histograms == TRUE){
        
            if (j == 1){
                h <- a$counts
            } else {
                h <- h + a$counts
            } 
     
            if (j %% 3 == 0){
                out <- h/3
        
                x <- c(out_matrix[j-2,],out_matrix[j-1,],out_matrix[j,]) 
 
                print(paste("DATA:",columns[i],j,mean(x),sd(x)))
                lines(a$mids,h/sum(h),lw=2,col=colors[j])
                h <- rep(0,length(h))
            } else {
                lines(a$mids,a$counts/sum(a$counts),lw=2,col=colors[i])
            }
        }

        }


    }

    # Create a legend
    names <- c()
    for (j in 1:length(files)){
        x <- strsplit(files[j],"_")
        y <- strsplit(x[[1]][3],"-")[[1]]

        names[j] <- paste(x[[1]][1],x[[1]][2],y)        

    }

    legend("topright",legend=names,col=colors,lw=2,lt=1,
           inset=0.05,bty='n')

    dev.off() 

}



