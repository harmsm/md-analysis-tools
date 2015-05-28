residues_to_take <- seq(1,250)
frame_cutoff <- 0

data_list <- c("SR2_N0600_run000", "SR2_N0600_run001")

do_barplot <- function(data,residue_list,frame_cutoff=800){

    d <- read.table(data)

    residues_to_take <- residue_list + 1
    frames_to_take <- d$frame > frame_cutoff

    d <- d[frames_to_take,residues_to_take]

    residues <- names(d)
    counts <- matrix(NA,4,length(residues))
    for (i in 1:length(residues)){
        v <- subset(d,select=residues[i])
        for (j in 1:4){
            counts[j,i] = length(v[v == j])
        }
    }

    counts <- counts/length(d[,1])

    # Grouped Bar Plot
    barplot(counts, main="",xlab="residue",
            col=c("red","orange","darkgreen","blue"),
            beside=FALSE,border=NA,space=0,axes=FALSE,ann=FALSE)
    #        legend=c("coil","turn","extended","helix"), beside=FALSE)

    counts

}

for (d in data_list){

    print(d)
    root <- d

    pdf(paste(root,".pdf",sep=""),height=3,width=10)

    x <- do_barplot(paste(d,"-processed.txt",sep=""),residues_to_take,
                    frame_cutoff)
    mtext(paste(d),3,line=1.5)
    dev.off()
}


