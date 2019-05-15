library(ggseqlogo)
argv <- commandArgs(TRUE)
print(argv)

f <- read.csv(as.character(argv[1]),header=FALSE,row.names=1)
x <- strsplit(as.character(argv[1]),"[.]")
node <- x[[1]][length(x[[1]])]
print(node)

svg(file = paste(list(node,"_logo.svg"),collapse = ""))
ggseqlogo(as.matrix(f),method="prob",col_scheme="clustalx")
dev.off()


    
