

######     HapNetComplexity.R        ######   ><(((º> 

# Haplotype Network Complexity Script #
# estimates complexity of haplotype networks 
# by combining indices of genetic and topological diversity

# Citation: 
# Garcia E, Wright D, Gatins R, Roberts MB, Pinheiro HT, Salas E, et al. (2021) Haplotype network branch diversity, a new metric combining genetic and topological diversity to compare the complexity of haplotype networks. PLoS ONE 16(6): e0251878. https://doi.org/10.1371/journal.pone.0251878
# Contact: Eric Garcia, e1garcia@odu.edu

# Use this script to build haplotype networks and calculate:

# Nucleotide diversity (Pi)                 = a measurement of the genetic distance between sequences
# Haplotype diversity (Hd)                  = a measurement of the genetic diversity in a population
# Branch diversity (Bd)                     = a new measurement of the topological diversity in a haplotype network or the diversity of interrelationships among the observed haplotypes in a population
# Haplotype Network Branch diversity (HBd)  = a new measurement of the complexity in a haplotype network

# Additional calculations:

# Number of individuals (n)
# Number of haplotypes (nH)
# Number of haplotype classes (nHc)
# Frequency of haplotype classes (niHc)

# See paper for definitions

# You will need:

# 1.- a fasta file with your sequences
# 2.- a file with site info (optional). 
    # The site file consists of a .csv file with the sample names in the first column and site names in the second. 
    # Call these columns "sample" and "site", respectively. Save the file as "sites.csv"
        # All metrics can still be calculated without any site information (see below) 



## First, create a haplotype network ##

# Install packages
install.packages('pegas')
install.packages('geiger')
install.packages('adegenet')
install.packages('reshape')

# Load packages
library(pegas)
library(geiger)
library(reshape)

# Clear your global environment if needed
rm(list=ls())

# If this script and data files live in the same directory
# Set the script's directory as the working directory  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Otherwise, specify the directory of your data files
#setwd("/Users/Yourname/YourDirectory")

# Load the alignment
#alignment<-read.dna("yourfile.fasta", format="fasta")

# Example
alignment <- read.dna("S16_File_Parupeneus_barberinus.fasta", format="fasta")
alignment <- read.dna("S25_File_testing_S1_D.fasta", format="fasta")

# Summary of the alignment
alignment

# Load or create site data frame:

# Load: 
# If you already have a .csv file for sites, load it as 
#sites<-read.csv('filename.csv')

#Example
#sites <- read.csv("S17_File_Parupeneus_barberinus_sites.csv")
#sites

# Create:
# Since our metrics are not affected by site information, 
# it is possible to create a data frame giving sites the same names 
# as the samples, making each site unique. 
# This guarantees that every sequence will always have a site match allowing the script to run.
# The following code will create this data frame. 
# Do not run this line if you have already generated a sites object above or it will be overwritten.

sites<- data.frame(sample=labels(alignment),site=labels(alignment))
sites

# The sample names should be an exact match in the fasta and site files.
# Check that this is true with this code: 
identical(as.character(sites$sample),labels(alignment))

# Calculate haplotype frequencies 
h <- haplotype(alignment)

#check haplotypes
h

# Build the haplotype network based on the haplotype frequencies
net <- haploNet(h)

# Summary of your network
net

# If you would like to see the list of the branches/links and relationships (step = no. of mutations) use print.default()
#print.default(net)

# Add column for haplotypes in sites
sites <- cbind(sites, haplotype=rep(NA, nrow(sites)))

# Assign haplotypes to the samples automatically using a loop
for (i in 1:length(labels(h))){
  sites$haplotype[attr(h, "index")[[i]]]<-i
}

# Build a table of frequencies that will tell R how to fill pies on haplotype network
ft <- table(sites[3:2])
ft

# Plot the haplotype network
plot(net, size=attr(net, "freq"), scale.ratio = 1, cex = 0.8, labels=TRUE, 
      pie = ft, font=2, fast=TRUE, threshold=0, show.mutation = 0, legend=FALSE)


# We have decided not to show the mutation steps between haplotypes to avoid clutter in our figures
# You can show mutations with a tick mark or a dot by setting "show.mutation" to 1 or 2, respectively
# Here, legend=FALSE since we give an improved customized legend option below 
# If you like to manually rearrange your plot, explore the following commands:
    # Be sure to hit Esc once you are done rearranging the nodes
    #P <- replot()
    #replot(P)


# Place a legend if desired

# Run the first command below, 
# then click on the plot area exactly where you would like to place the legend, 
# and finally, run the second command to draw it

# Set the coordinates for legend
#legend_coord <- locator(1)
#plot the legend
#legend(x=legend_coord[1], y=legend_coord[2], unique(sort(sites$site)),
#       fill=rainbow(length(unique(sites$site))), bty="n", cex=1)



### Now, we can calculate diversity metrics ###


## Calculate Nucleotide diversity (Pi) ##     ><(((º> 

Pi <- nuc.div(alignment)
Pi


## Calculating Haplotype diversity (Hd) ##    

# Transform you frequency table into a data frame to calculate haplotype statistics 
ft<- as.data.frame.matrix(ft)


# Add total column to ft and sum the rows
ft$total <- apply(ft, 1, sum)


# Add a column and calculate the haplotype square frequencies 
ft<-cbind(ft, sqfreq = (ft$total/sum(ft$total))^2) 
ft

# Sum the square frequencies
Efh2 <- sum(ft$sqfreq)

# Retrieve total no. of individuals
n <- nrow(sites)
n

## Calculate Haplotype diversity (Hd) ##      ><(((º>

Hd <- (1-Efh2)*(n/(n-1))
Hd


## Calculating Branch diversity (Bd) ##      

# Retrieve all haplotype connections
haplotype_counts <-c(net[,1],net[,2])

# Make a data frame listing the haplotypes in the first column and their corresponding no. of branches in the second
Hbranches <- aggregate(data.frame(branches = haplotype_counts), list(haplotype = haplotype_counts), length)
Hbranches

# Add a column to note the number of individuals per haplotype
Hbranches <- cbind(Hbranches, nseq = ft$total)
Hbranches

# Count the total number of haplotypes (nH)
nH <- nrow(Hbranches)
nH

# Aggregate by Haplotype class (Hc) (classes are set according to their number of branches)
Hclasses_temp <- aggregate(Hbranches,list("Hc" = Hbranches$branches), sum)

# Keep only the list of haplotype classes (Hc) and their frequencies (i.e. number of individual, niHc). Ignore the middle columns.
Hclasses <- Hclasses_temp[,c(1,4)]    
names(Hclasses)[2] <- "niHc"
Hclasses  # Hc and niHc

# Count the total number of haplotype classes (nHc)
nHc <- nrow(Hclasses)
nHc

# Calculate class square frequencies (fHc2)
Hclasses<-cbind(Hclasses, sqfreq = (Hclasses$`niHc`/n)^2) 

# Sum of class square frequencies (∑fnHi2)
EfHc2 <- sum(Hclasses$sqfreq)

## Calculate Branch diversity (Bd) ##     ><(((º>
Bd <- (1-EfHc2)*(n/(n-1))
Bd

# Calculate Haplotype network branch diversity (HBd)    ><(((º>   --¥-¥-{

HBd<-(1-Efh2)*(1-EfHc2)*(n/(n-1))
HBd
    
