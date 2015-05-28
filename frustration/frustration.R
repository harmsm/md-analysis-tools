# Has a zero-padding string function (pad0)
library(PBSmodelling)

read_data <- function(root="SR2_N0600",type="dihed-raw",frame_cutoff=0,
                      run_list=c(0,1,2)){

    # Function to read in data of a given type for a given set of replicate
    # runs.  The function also tosses frames from the burn in.

    f1 <- paste("trajectory-data/",root,"_run",sep="")
    f3 <- paste("_",type,"-processed.txt",sep="")

    for (i in 1:length(run_list)){
        tmp <- read.table(paste(f1,pad0(run_list[i],3),f3,sep=""))
        tmp$run <- rep(run_list[i],length(tmp[,1]))
        if (i == 1){
            d <- tmp 
        } else {
            d <- rbind(d,tmp)
        }
    }

    d <- subset(d,frame>frame_cutoff)

    id <- subset(d,select=c(run,frame))
    not_id <- subset(d,select=-c(run,frame))

    cbind(id,not_id)

}

get_r41_rotamer <- function(chi_table,class="gln"){

    # Place r41 into roamter bins based on the chi1, chi2, and chi3 values 
    # specified in the "bin" list.  This is entirely empirical, designed based 
    # on the clusters of populated chi1/chi2/chi3 values over all trajectories.

    if (class == "gln"){

        bin_names <- c("a","b","c","d","e","f","ga","gb","h","i")
   
        bin <- list() 
        bin[[1]]  <- c(0  ,120,120,240,120,360)
        bin[[2]]  <- c(120,240,  0,120,  0,180)

        bin[[3]]  <- c(240,360,120,240,180,300)

        bin[[4]]  <- c(120,240,120,240,  0,180)
        bin[[5]]  <- c(120,240,  0,180,180,360)

        bin[[6]]  <- c(240,360,120,240, 60,180)
        bin[[7]]  <- c(240,360,120,240,  0, 60)
        bin[[8]]  <- c(240,360,120,240,300,360)
    
        bin[[9]]  <- c(240,360,240,360,180,360)
        bin[[10]] <- c(240,360,240,360,  0,180)

    } else if (class == "glu"){

        bin_names <- c("v","w","x","y","z")

        bin <- list()
        bin[[1]] <- c(  0,120,120,240,  0,360)
        bin[[2]] <- c(120,240,  0,120,  0,360)
        bin[[3]] <- c(120,240,120,240,  0,360)
        bin[[4]] <- c(240,360,120,240,  0,360)
        bin[[5]] <- c(240,360,240,360,  0,360)

    }

    bin_output <- rep(NA,length(chi_table[,1]))
    for (i in 1:length(chi_table[,1])){
        print(i)
        for (j in 1:length(bin)){
            if ((chi_table[i,1] >= bin[[j]][1] & chi_table[i,1] < bin[[j]][2]) &
                (chi_table[i,2] >= bin[[j]][3] & chi_table[i,2] < bin[[j]][4]) &
                (chi_table[i,3] >= bin[[j]][5] & chi_table[i,3] < bin[[j]][6])){
                bin_output[i] <- bin_names[j] 
                break
            }
        }
    }
  
 
    bin_output

}

create_column_subset <- function(d,column_list){

    sel <- intersect(names(d),column_list)
    apply(subset(d,select=sel),1,sum)

}

assign_micro_bins <- function(root="SR2_N0600",class="gln",
                              frame_cutoff=0,run_list=c(0,1,2)){

    # This function reads in the relevant information for the "root" replicate
    # trajectories, then places each frame into a micro_bin by combining the 
    # Gln41 rotamer and a set of recorded hydrogen bonds.

    # Read in dihedral data and classify based on r41_rotamer
    d <- read_data(root,type="dihed-raw",frame_cutoff,run_list)
    r41 <- subset(d,select=c(r41_chi1,r41_chi2,r41_chi3))
    d$r41_rotamer <- get_r41_rotamer(r41,class)

    # Read in hydrogen bond data
    e <- read_data(root,type="hbond",frame_cutoff,run_list)
    e <- subset(e,select=-c(run,frame))
    
    # Read in internal water data
    f <- read_data(root,type="water",frame_cutoff,run_list)
    f <- subset(f,select=-c(run,frame))

    # Put everyone together
    d <- cbind(d,e,f)
    
    # Toss unclassified rotamers
    d <- d[!is.na(d$r41_rotamer),]
    
    if (class == "gln"){        
    
        group_1 <- c("a")
        group_2 <- c("b")
        group_3 <- c("c")
        group_4 <- c("d")
        group_5 <- c("e")
        group_6 <- c("f")
        group_7 <- c("ga","gb")
        group_8 <- c("h")
        group_9 <- c("i")
        
        rotamer_groups <- list(group_1,group_2,group_3,group_4,group_5,group_6,
                               group_7,group_8,group_9)

    } else if (class == "glu") {
        
        group_1 <- c("v")
        group_2 <- c("w")
        group_3 <- c("x")
        group_4 <- c("y")
        group_5 <- c("z")

        rotamer_groups <- list(group_1,group_2,group_3,group_4,group_5)

    } else {

        print(paste("Class",class,"not recognized!"))
        return(NA)

    }
    
    rotamer <- rep(NA,length(d[,1]))
    for (i in 1:length(d[,1])){
        for (j in 1:length(rotamer_groups)){
            if (is.element(d$r41_rotamer[i],rotamer_groups[[j]])){
                rotamer[i] <- j
                next
            }
        }
    }


    # Some important hydrogen bonds to count
    ne2_hoh <- create_column_subset(d,c("NE2.41GLN_SOL.ante",
                                        "NE2.41GLN_SOL.lower",
                                        "NE2.41GLN_SOL.upper",
                                        "NE2.41GLN_SOL.inter",
                                        "NE2.41GLN_SOL.bulk"))

    ne2_lig <- create_column_subset(d,c("NE2.41GLN_O3.249LIG"))
   
    hoh_oe <- create_column_subset(d,c("SOL.ante_OE.41GLU",
                                       "SOL.lower_OE.41GLU",
                                       "SOL.upper_OE.41GLU",
                                       "SOL.inter_OE.41GLU",
                                       "SOL.bulk_OE.41GLU",
                                       "SOL.ante_OE1.41GLN",
                                       "SOL.lower_OE1.41GLN",
                                       "SOL.upper_OE1.41GLN",
                                       "SOL.inter_OE1.41GLN",
                                       "SOL.bulk_OE1.41GLN"))
    
    lig_oe   <- create_column_subset(d,c("O3.249LIG_OE1.41GLN",
                                         "O3.249LIG_OE.41GLU"))
    
    hoh_lig <- create_column_subset(d,c("SOL.ante_O3.249LIG",
                                        "SOL.lower_O3.249LIG",
                                        "SOL.upper_O3.249LIG",
                                        "SOL.inter_O3.249LIG",
                                        "SOL.bulk_O3.249LIG"))

    lig_hoh <- create_column_subset(d,c("O3.249LIG_SOL.ante",
                                        "O3.249LIG_SOL.lower",
                                        "O3.249LIG_SOL.upper",
                                        "O3.249LIG_SOL.inter",
                                        "O3.249LIG_SOL.bulk"))
    
    r82_r41 <-  create_column_subset(d,c("NH.82ARG_OE.41GLU",
                                         "NH.82ARG_OE1.41GLN"))

    r41_r75 <- create_column_subset(d,c("NE2.41GLN_SD.75MET"))

    d$micro_bin <- paste(rotamer,
                         ne2_hoh,
                         ne2_lig,
                         hoh_oe, 
                         lig_oe,
                         r41_r75,
                         r82_r41,
                         sep="-")

    d

}

meta_reader <- function(){

    n06 <- assign_micro_bins(root="SR2_N0600",class="gln")
    el_n06 <- assign_micro_bins(root="SR2-Q41E-M75L_N0600",class="glu")
    n15 <- assign_micro_bins(root="SR2_N1500",class="gln")
    el_n15 <- assign_micro_bins(root="SR2-Q41E-M75L_N1500",class="glu")
    npr <- assign_micro_bins(root="SR2_NPR",class="gln")

    d <- list()
    d$n06 <- n06
    d$el_n06 <- el_n06
    d$n15 <- n15
    d$el_n15 <- el_n15
    d$npr <- npr

    d
}


get_trans_matrix <- function(d,condition_on_population=FALSE){

    # Determine a transition matirx between micro bins

    # Convert the string micro bin names into a vector of integers
    all_bins <- levels(factor(d$micro_bin))
    num_bins <- length(all_bins)

    bin_dict <- list()
    for (i in 1:num_bins){
        bin_dict[[all_bins[i]]] <- i
    }

    micro_bin <- unlist(bin_dict[d$micro_bin])

    # Walk through each step, recording each transition.  Make sure that the
    # frames are contiguous; otherwise, skip the step.
    m <- matrix(0,num_bins,num_bins)
    old_bin <- micro_bin[1]
    for (i in 2:length(micro_bin)){
    
        if (d$frame[i] == (d$frame[i-1] + 1)){
            
            m[old_bin,micro_bin[i]] = m[old_bin,micro_bin[i]] + 1

            if (old_bin != micro_bin[i]){
                old_bin = micro_bin[i]
            }
        } else {
            old_bin <- micro_bin[i]
        }

    }
    
    # Divide each row by the population of that micro bin
    if (condition_on_population == TRUE){
        for (i in 1:num_bins){
            m[i,] <- m[i,]/apply(m,2,sum)
        }
    }
        

    m

}

write_micro_bins <- function(d,root="SR2_N0600"){

    # Write out the micro_bin of every frame in the input trajectories.

    bins <- levels(factor(d$micro_bin))

    for (s in bins){
    
        to_write <- subset(d,micro_bin==s)
        to_write <- cbind(to_write$run,to_write$frame,to_write$micro_bin)

        sink(paste(root,"_bin_",s,".txt",sep=""))
        print(to_write)
        sink()

    }
}



hbond_stats <- function(d){

    source('useful-selections.R')
    
    micro_bins <- levels(factor(d$micro_bin))

    out <- as.data.frame(rep(0,length(micro_bins)))

    names(out) <- c("micro_bin")
    out$r75_acc_m <- rep(0,length(micro_bins))
    out$r75_acc_s <- rep(0,length(micro_bins))
    out$r41_acc_m <- rep(0,length(micro_bins))
    out$r41_acc_s <- rep(0,length(micro_bins))
    out$r41_don_m <- rep(0,length(micro_bins))
    out$r41_don_s <- rep(0,length(micro_bins))
    out$lig_acc_m <- rep(0,length(micro_bins))
    out$lig_acc_s <- rep(0,length(micro_bins))
    out$lig_don_m <- rep(0,length(micro_bins))
    out$lig_don_s <- rep(0,length(micro_bins))
    out$r82_don_m <- rep(0,length(micro_bins))
    out$r82_don_s <- rep(0,length(micro_bins))


    for (i in 1:length(micro_bins)){
        out$micro_bin[i] <- micro_bins[i]

        a <- subset(d,micro_bin==micro_bins[i])
    
        out$r75_acc_m[i] <- mean(apply(subset(a,select=r75_acc),1,sum))
        out$r75_acc_s[i] <- sd(apply(subset(a,select=r75_acc),1,sum))
        
        out$r41_acc_m[i] <- mean(apply(subset(a,select=r41_acc),1,sum))
        out$r41_acc_s[i] <- sd(apply(subset(a,select=r41_acc),1,sum))
        out$r41_don_m[i] <- mean(apply(subset(a,select=r41_don),1,sum))
        out$r41_don_s[i] <- sd(apply(subset(a,select=r41_don),1,sum))
        
        out$lig_acc_m[i] <- mean(apply(subset(a,select=lig_acc),1,sum))
        out$lig_acc_s[i] <- sd(apply(subset(a,select=lig_acc),1,sum))
        out$lig_don_m[i] <- mean(apply(subset(a,select=lig_don),1,sum))
        out$lig_don_s[i] <- sd(apply(subset(a,select=lig_don),1,sum))

        out$r82_don_m[i] <- mean(apply(subset(a,select=r82_don),1,sum))
        out$r82_don_s[i] <- sd(apply(subset(a,select=r82_don),1,sum))

    }

    out
}


micro_bin_summaries <- function(d,to_take=c(),cutoff=NA,remove_zero=FALSE){

    micro_bins <- levels(factor(d$micro_bin))

    if (length(to_take) == 0){
        to_skip <- c("macro_bin","run","frame","rotamer",
                     "r41_rotamer",
                     names(d)[grep("chi",names(d))])
        to_take <- setdiff(names(d),to_skip)
    }

    a <- d[,to_take]

    counts <- c()
    taken <- c()
    out <- matrix(NA,length(micro_bins),length(a[1,])-1)
    for (i in 1:length(micro_bins)){

        b <- subset(a,micro_bin==micro_bins[i],select=-micro_bin)

        if (! is.na(cutoff)){
            if (length(b[,1]) < cutoff){
                counts[i] <- NA
                next
            }
        }

        taken <- c(taken,i)
        out[i,] <- apply(b,2,mean)
        counts[i] <- length(b[,1]) 
    }

    out <- as.data.frame(out)
    names(out) <- names(a)[1:length(a)-1]
    
    if (remove_zero){
        col <- apply(out,2,sum) != 0
        out <- out[,col]

    }

    out$micro_bin <- micro_bins
    out$counts <- counts

    out <- out[taken,] 


    out


}

transition_summaries <- function(d,cutoff=-1){

    m <- get_trans_matrix(d)
    micro_bin_names <- levels(factor(d$micro_bin))

    micro_bin_totals <- apply(m,2,sum)
    ord <- order(micro_bin_totals,decreasing=TRUE)

    print("State populations:")
    pop <- cbind(ord,micro_bin_names[ord],micro_bin_totals[ord],
                 micro_bin_totals[ord]/sum(micro_bin_totals))

    print(pop)    

    for (i in 1:length(pop[,1])){

        if (pop[i,3] > cutoff){
            
            index <- as.numeric(pop[i,1])

            trans <- m[index,]
            ord <- order(m[index,],decreasing=TRUE)               
 
            print(paste(pop[i,2],pop[i,4]))
            out <- cbind(ord,micro_bin_names[ord],trans[ord],
                        trans[ord]/(micro_bin_totals[index]))

            print(out[out[,3] != 0,])

             
        }
    }

    pop

}


analyze_macro_bins <- function(d,macro_bin_list){

    to_skip <- c("macro_bin","micro_bin","run","frame","rotamer",
                 "r41_rotamer")
    to_skip <- c(to_skip,names(d)[grep("chi",names(d))])

    d$macro_bin <- unlist(macro_bin_list[d$micro_bin])

    all_bins <- levels(factor(d$macro_bin))
    num_bins <- length(all_bins)

    bin_dict <- list()
    for (i in 1:num_bins){
        bin_dict[[as.numeric(all_bins[i])]] <- i
    }

    d$macro_bin <- unlist(bin_dict[d$macro_bin])

    m <- matrix(0,num_bins,num_bins)

    old_micro_bin <- d$micro_bin[1]
    old_macro_bin <- d$macro_bin[1]
    for (i in 2:length(d[,1])){
  
        # Make sure the frames are contiguous; otherwise, skip the step
        if (d$frame[i] == (d$frame[i-1] + 1)){
            
            m[old_macro_bin,d$macro_bin[i]] = m[old_macro_bin,d$macro_bin[i]] + 1

            if (old_macro_bin != d$macro_bin[i]){
                old_macro_bin = d$macro_bin[i]
                old_micro_bin = d$micro_bin[i]
            }

        } else {
            old_macro_bin <- d$macro_bin[i]
            old_micro_bin = d$micro_bin[i]
        }

    }

    state_pops <- apply(m,1,sum)
    barriers <- -log(m/state_pops)
    dG <- -log(state_pops/state_pops[1])
  
    output <- list()
    output[["trans_mat"]] <- m
    output[["state_pops"]] <- state_pops
    output[["barriers"]] <- barriers
    output[["dG"]] <- dG
    output[["bins"]] <- all_bins


    output
}


macro_bin_interactions <- function(d,macro_bin_list,cutoff=-1){
    
    to_skip <- c("macro_bin","micro_bin","run","frame","rotamer",
                 "r41_rotamer")
    to_skip <- c(to_skip,names(d)[grep("chi",names(d))])

    d$macro_bin <- unlist(macro_bin_list[d$micro_bin])
    all_bins <- levels(factor(d$macro_bin))
    print(all_bins)

    to_take <- setdiff(names(d),to_skip)
    interactions <- list()
    for (b in all_bins){
        a <- subset(d,macro_bin==b,select=to_take)
        x <- apply(a,2,mean)
        interactions[[b]] <- x[x>cutoff]
    }

    to_check <- c("NH.82ARG_","_SD.75MET","_OE1.41GLN","NE2.41GLN_","_OE.41GLU",
                  "O3.249LIG_","_O3.249LIG","SOL.bulk","SOL.ante","SOL.lower",
                  "SOL.inter","SOL.upper")

    for (b in all_bins){
        for (c in to_check){
            print(paste(b,c,sum(interactions[[b]][grep(c,names(interactions[[b]]))])))
        }
    }
        

    interactions
}


qual <- function(m,inside){

    outside <- setdiff(seq(1,length(m[1,])),inside)

    # Take transitions within this matrix
    within_macrobin_matrix <- m[inside,inside]

    # Self transitions are along the diagnol
    to_self <- sum(diag(within_macrobin_matrix))

    # To others within the macrobin are off diagnol
    to_inside <- sum(within_macrobin_matrix) #- to_self

    # All other transfers
    outside_macrobin_matrix <- m[inside,outside]
    
    to_outside <- sum(outside_macrobin_matrix)

    to_inside/to_outside

}

cluster_quality <- function(d,macro_bin_list,relative_counts=FALSE){

    m <- get_trans_matrix(d)
    if (relative_counts){
        m <- m/apply(m,2,sum)
    }

    x <- macro_bin_list[levels(factor(d$micro_bin))]
    y <- as.data.frame(unlist(x))

    names(y) <- c("macro_bin")
    macro_bins <- levels(factor(y$macro_bin))
    micro_bins <- levels(factor(d$micro_bin))

    clusters_out <- list()
    out <- matrix(NA,length(macro_bins),3)
    for (i in 1:length(macro_bins)) { 
       
        c <- seq(1:length(m[1,]))[y$macro_bin == macro_bins[i]]

        out[i,1] <- qual(m,c)
        out[i,2] <- sum(m[,c])/sum(m)
        out[i,3] <- length(c)

        clusters_out[[i]] <- c

        
    }

    o <- order(out[,2],decreasing=TRUE)
    out <- out[o,]
    macro_bins <- macro_bins[o]
    clusters_out <- clusters_out[o]
    num_frames <- length(d$micro_bin)
    cum_sum <- cumsum(out[,2]) 

    for (i in 1:length(macro_bins)){

        print(sprintf("Cluster %s: (%i) %8.3f %8.3f %8.3f %8.0f %8i",
                      macro_bins[i],i,out[i,1],out[i,2],cum_sum[i],
                      (out[i,2]*num_frames),out[i,3]))
        print(micro_bins[clusters_out[[i]]])

    }

}


write_total_summary <- function(d,macro_bin_list,outfile_root){

    pdf(paste(outfile_root,".pdf",sep=""),height=8,width=10)
    o <- plot_cluster(d,trajectory=outfile_root)
    text(seq(1,length(macro_bin_list)),y=-1,
         labels=unlist(macro_bin_list[o]),col='red')
    dev.off()
   
    sink(paste(outfile_root,".txt",sep=""))
 
    print("*************************************************************")
    print("Cluster quality:")
    print(cluster_quality(d,macro_bin_list,raw_counts=TRUE))

    print("*************************************************************")
    print("Cluster summary:")
    print(analyze_macro_bins(d,macro_bin_list))

    print("*************************************************************")
    print("Cluster interactions:")
    print(macro_bin_interactions(d,macro_bin_list,cutoff=0))

    print("*************************************************************")
    print("Micro bin transitions:") 
    transition_summaries(d)

    print("*************************************************************")
    print("Micro bin interactions:")
    print(micro_bin_summaries(d,remove_zero=TRUE))

    sink()

}



tarjan_core <- function(m,threshold=0){

    # Find strongly connected components of a graph using the Tarjan 
    # algorithm.  
    #   m: transition matrix (in dG; larger is *worse*)
    #   threshold: maximum free energy to score nodes as connected.

    # A set of variables used to track whether a given vertex is already in a
    # component, etc.
    d <- list()
    d$m <- m
    d$threshold <- threshold
    d$num_vert <- length(m[1,])
    d$index <- 0
    d$stack <- c()
    d$vert_index <- rep(NA,d$num_vert)
    d$vert_lowlink <- rep(NA,d$num_vert)
    d$SSC_counter <- 1
    d$SSCs <- list()

    strongconnect <- function(i,d){

        # Recursive function to do a depth-first search of vertices connected to
        # vertex i
  
        d$vert_index[i] <- d$index
        d$vert_lowlink[i] <- d$index 
        d$index <- d[["index"]] + 1

        # Push i onto the stack
        d$stack <- c(d$stack,i)

        for (j in 1:d$num_vert){

            # Don't make self-self comparison
            if (j == i) { 
                next 
            }

            # Only consider components connected by value higher than threshold
            if (d$m[i,j] >= d$threshold){
                next
            }

            if (is.na(d$vert_index[j])){
                # If vertex j has not been seen, call strong connect
                d <- strongconnect(j,d)
                d$vert_lowlink[i] <- min(c(d$vert_lowlink[i],d$vert_lowlink[j]))
            } else if (is.element(j,d$stack)){
                # Succesor j is in the stack, thus in the current SSC
                d$vert_lowlink[i] <- min(c(d$vert_lowlink[i],d$vert_index[j]))
            }
        }

        # i is a root node, we've found a complete SSC; write it out to SSCs
        if (d$vert_lowlink[i] == d$vert_index[i]){

            # Pop the stack
            j = d$stack[length(d$stack)]
            if (length(d$stack) == 1){
                d$stack <- c()
            } else {
                d$stack <- d$stack[1:(length(d$stack)-1)]
            }
            
            # Start a new strongly connected component
            d$SSCs[d$SSC_counter] <- c(j)
   
            while (i != j){
 
                # Pop the stack
                j = d$stack[length(d$stack)]
                if (length(d$stack) == 1){
                    d$stack <- c()
                } else { 
                    d$stack <- d$stack[1:(length(d$stack)-1)]
                }

                # Append j to the current SSC
                d$SSCs[[d$SSC_counter]] <- c(d$SSCs[[d$SSC_counter]],j)

            }
    
            # Increment counter
            d$SSCs[[d$SSC_counter]] <- unique(d$SSCs[[d$SSC_counter]])
            d$SSC_counter <- d$SSC_counter + 1 
              
        }

        return(d)

    }

    # Run through all vertices, starting a strongconnect recursion if the 
    # vertex has not been visited yet.
    for (i in 1:d$num_vert){
        if (is.na(d$vert_index[i])){
            d <- strongconnect(i,d)    
        } 
    }

    d$SSCs

}


dfs_core <- function(m,threshold=0){

    # Find isolated subgraphs using depth-first search
    #   m: transition matrix (in dG; larger is *worse*)
    #   threshold: maximum free energy to score nodes as connected.

    # A set of variables used to track whether a given vertex is already in a
    # component, etc.
    d <- list()
    d$m <- m
    d$threshold <- threshold
    d$num_vert <- length(m[1,])
    d$vert_seen <- rep(FALSE,d$num_vert)
    d$cluster_counter <- 1
    d$clusters <- list()

    connect <- function(i,d){

        # Recursive function to do a depth-first search of vertices connected to
        # vertex i

        for (j in 1:d$num_vert){

            # Don't make self-self comparison
            if (j == i) { 
                next 
            }

            # Only consider components connected by value higher than threshold
            if (d$m[i,j] >= d$threshold){
                next
            }

            if (!d$vert_seen[j]){
                # If vertex j has not been seen, record that it is part of this
                # cluster and then look at its leaves
                d$vert_seen[j] <- TRUE
                d$clusters[[d$cluster_counter]] <- c(d$clusters[[d$cluster_counter]],j)
                d <- connect(j,d)
            }
        }

        return(d)

    }

    # Run through all vertices, starting a connect recursion if the 
    # vertex has not been visited yet.
    for (i in 1:d$num_vert){
        if (! d$vert_seen[i]){
            d$vert_seen[i] <- TRUE
            d$clusters[[d$cluster_counter]] <- c(i)
            d <- connect(i,d)
            d$cluster_counter <- d$cluster_counter + 1
        } 
    }

    d$clusters

}

subgraph_core <- function(m,threshold=0){
 
    # Cluster the matrix using subgraphs

    size <- length(m[1,]) 

    # Initialize each microbin as its own cluster
    micro_to_macro <- list()
    for (i  in 1:size){
        micro_to_macro[i] <- i
    }

    # Go through the distance matrix and use the threshold to decide whether
    # two states belong together or not.
    seen_clusters <- c()
    for (i in 1:size){
        for (j in 1:size){
            if (i == j) { next}
            if (m[i,j] < threshold){
                
                # Use seen_clusters to decide whether to assign i->j or j->i
                # so we don't end up with synonymous cluster names.
                if (is.element(micro_to_macro[[i]],seen_clusters)){
                    micro_to_macro[[j]] <- micro_to_macro[[i]]
                } else if (is.element(micro_to_macro[[j]],seen_clusters)){
                    micro_to_macro[[i]] <- micro_to_macro[[j]]
                } else {
                    micro_to_macro[[i]] <- micro_to_macro[[j]]
                    seen_clusters <- c(seen_clusters,micro_to_macro[[j]])
                }
            }
        }
    }

    # CURRENTLY BROKEN
    # Convert the state->cluster data structure to cluster->{state1,state2...}
    # structure
    #clusters <- list()
    #for (n in unique(unlist(micro_to_macro))){
    #    clusters[[paste(n)]] <- seq(1,length(micro_to_macro)[unlist(micro_to_macro) == n])
    #}
    #print(clusters)  

    NA
}

graph_cluster <- function(d,method="tarjan",threshold=0,collapse=NA,T=297){

    # d is a set of data for a ligand/protein pair
    # threshold is the free energy threshold below which to consider two nodes
    #   connected.
    # method is one of the following:
    #   "tarjan": Use the Tarjan algorithm to find strongly connected components
    #   "dfs": Use a simple depth-first search to find connected nodes.
    #   "subgraph": Use a simple subgraph cutting algorithm to find connected 
    #       nodes.  
    # collapse: a function applied to decide whether to take (i,j) or (j,i) if
    #   you want to convert an asymmetrical to a symmetrical matrix.  "subgraph"
    #   requires a collapse function.  If none is specified, it uses "min."
    # T: temperature for converting a the transition probability from "d" into
    #   a free energy.

    # Figure out what clustering method to use
    methods <- list()
    methods[["tarjan"]] <- tarjan_core
    methods[["subgraph"]] <- subgraph_core
    methods[["dfs"]] <- dfs_core
    
    method_fcn <- try(methods[[method]])
    if (!is.function(method_fcn)){
        print("You must specify one of the following clustering methods:")
        print(names(methods))
        return(NA)
    }  

    # Define a collpase method for subgraph if this was not specified
    if (is.na(collapse) & method == "subgraph"){
        print("No collapse function specified for subgraph.  Using min.")
        collapse <- min
    }

    # Create transition matrix and grab names of the populated micro bins. 
    raw_trans_matrix <- get_trans_matrix(d)
    micro_bins <- levels(factor(d$micro_bin))

    # Order the matrix from largest to smallest value
    o <- order(apply(raw_trans_matrix,2,sum),decreasing=TRUE)
    raw_trans_matrix <- raw_trans_matrix[o,o]
    micro_bins <- micro_bins[o]
    
    # Convert these transitions to free energies
    trans_matrix <- raw_trans_matrix/apply(raw_trans_matrix,1,sum)
    trans_matrix <- -(0.001987)*T*log(trans_matrix)
  
    # Convert any NA values to infinite free energies 
    diag(trans_matrix) <- 0
    size <- length(trans_matrix[1,])
    for (i in 1:size){
        for (j in 1:size){
            if (is.na(trans_matrix[i,j])) { trans_matrix[i,j] <- Inf } 
        }
    }

    # Generate clusters
    clusters <- method_fcn(trans_matrix,threshold)

    # Create an output dictionary keying micro_bins to clusters
    out_list <- list()
    repeated <- c()
    for (i in 1:length(clusters)){
        for (j in 1:length(clusters[[i]])){
            repeated <- c(repeated,micro_bins[clusters[[i]][[j]]])
            out_list[[micro_bins[clusters[[i]][j]]]] <- i
        }
    }

    out_list 

}

cluster_cumplot <- function(data,method="tarjan",threshold=0,xlim=c(0,10)){


    print("el/N06")
    a <- analyze_macro_bins(data$el_n06,
                            graph_cluster(data$el_n06,method,threshold))
    print("QM/N06")
    b <- analyze_macro_bins(data$n06,
                            graph_cluster(data$n06,method,threshold))
    print("el/N15")
    c <- analyze_macro_bins(data$el_n15,
                            graph_cluster(data$el_n15,method,threshold))
    print("QM/N15")
    d <- analyze_macro_bins(data$n15,
                            graph_cluster(data$n15,method,threshold))
    print("QM/NPR")
    e <- analyze_macro_bins(data$npr,
                            graph_cluster(data$npr,method,threshold))


    a1 <- exp(-a$dG)/sum(exp(-a$dG))
    b1 <- exp(-b$dG)/sum(exp(-b$dG))
    c1 <- exp(-c$dG)/sum(exp(-c$dG))
    d1 <- exp(-d$dG)/sum(exp(-d$dG))
    e1 <- exp(-e$dG)/sum(exp(-e$dG))

    plot(NA,xlim=xlim,ylim=c(0,1),ann=FALSE,axes=FALSE)
    axis(1); mtext("macrostate, reverse rank ordered by size",1,line=2.5)
    axis(2); mtext("cumulative fraction of frames",2,line=2.5)

    lines(cumsum(a1[order(a1,decreasing=TRUE)]),col='orange',lw=2)
    lines(cumsum(b1[order(b1,decreasing=TRUE)]),col='red',lw=2)
    lines(cumsum(c1[order(c1,decreasing=TRUE)]),col='lightblue',lw=2)
    lines(cumsum(d1[order(d1,decreasing=TRUE)]),col='blue',lw=2)
    lines(cumsum(e1[order(e1,decreasing=TRUE)]),col='darkgreen',lw=2)
    
    points(cumsum(a1[order(a1,decreasing=TRUE)]),col='orange',pch=20,cex=1.2)
    points(cumsum(b1[order(b1,decreasing=TRUE)]),col='red',pch=20,cex=1.2)
    points(cumsum(c1[order(c1,decreasing=TRUE)]),col='lightblue',pch=20,cex=1.2)
    points(cumsum(d1[order(d1,decreasing=TRUE)]),col='blue',pch=20,cex=1.2)
    points(cumsum(e1[order(e1,decreasing=TRUE)]),col='darkgreen',pch=20,cex=1.2)

    legend("bottomright",c("el/N06","QM/N06","el/N15","QM/N15","QM/NPR"),
           lw=2,lt=1,col=c("orange","red","lightblue","blue","darkgreen"),
           inset=0.05,bty="n")


}


state_distrib <- function(data,method="tarjan",cutoff=0.95){

    threshold_list <- seq(0.0,3,by=0.1)
    labels <- c("el/NPT","QM/NPT","el/norP","QM/norP","QM/NPR")
    colors <- c("orange","red","lightblue","blue","darkgreen")


    out <- matrix(NA,length(threshold_list),5)
    for (i in 1:length(threshold_list)){

        print(threshold_list[i])
    
        a <- analyze_macro_bins(data$el_n06,
                                graph_cluster(data$el_n06,method,
                                threshold_list[i]))
        b <- analyze_macro_bins(data$n06,
                                graph_cluster(data$n06,method,
                                threshold_list[i]))
        c <- analyze_macro_bins(data$el_n15,
                                graph_cluster(data$el_n15,method,
                                threshold_list[i]))
        d <- analyze_macro_bins(data$n15,
                                graph_cluster(data$n15,method,
                                threshold_list[i]))
        e <- analyze_macro_bins(data$npr,
                                graph_cluster(data$npr,method,
                                threshold_list[i]))

        a1 <- exp(-a$dG)/sum(exp(-a$dG))
        b1 <- exp(-b$dG)/sum(exp(-b$dG))
        c1 <- exp(-c$dG)/sum(exp(-c$dG))
        d1 <- exp(-d$dG)/sum(exp(-d$dG))
        e1 <- exp(-e$dG)/sum(exp(-e$dG))

        out[i,1] <- which(cumsum(a1[order(a1,decreasing=TRUE)]) > cutoff)[1]
        out[i,2] <- which(cumsum(b1[order(b1,decreasing=TRUE)]) > cutoff)[1]
        out[i,3] <- which(cumsum(c1[order(c1,decreasing=TRUE)]) > cutoff)[1]
        out[i,4] <- which(cumsum(d1[order(d1,decreasing=TRUE)]) > cutoff)[1]
        out[i,5] <- which(cumsum(e1[order(e1,decreasing=TRUE)]) > cutoff)[1]

    }

    plot(NA,xlim=c(0,3),ylim=c(0,90),ann=FALSE,axes=FALSE)
    axis(1); mtext("energy cutoff (kcal/mol)",1,line=2.5)
    axis(2); mtext("number of states in top 95%",2,line=2.5)
   
    for (i in 1:5){
        points(threshold_list,out[,i],pch=20,cex=1.2,col=colors[i])
        lines(threshold_list,out[,i],lw=2,col=colors[i])
    }
  
    legend("topright",labels,lw=2,lt=1,col=colors,inset=0.05,bty="n")

    out

}


