
library(igraph)
library(lsa)
library(plyr)
library(compare)

#Initializing the data
graph_path <- "~/Documents/06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection/data/fb_caltech_small_edgelist.txt"
attr_path <- "~/Documents/06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection/data/fb_caltech_small_attrlist.csv"
attr <- read.csv(attr_path, header = TRUE)
graph <- read.graph(graph_path, format = "edgelist", directed = FALSE)
alpha <- 0.5

ans <- list()
for(i in 1:vcount(graph))
{
  ans[[i]] <- i
}

#Looping Phase 1 till convergence
repeat
{
  n <- vcount(graph)
  e <- ecount(graph)
  print(paste(n,e))
  community <- 1:n
  #Precompute the similarity matrix
  similarity <- matrix(ncol = n, nrow = n)
  for(i in 1:n)
  {
   for(j in i:n)
   {
     similarity[i,j] <- cosine(as.numeric(attr[i,]), as.numeric(attr[j,]))
     similarity[j,i] <- similarity[i,j]
   }
  }

  #Store the community before working on it
  prev <- community
  
  #Apply Phase 1
  community <- phase1(community, alpha)

  #Mapping the old vertices to new vertices..New Numbering
  mapped.communities <- mapvalues(community, from = sort(unique(community)), to = rep(1:length(unique(community)), each=1))
  
  graph.contracted <- contract.vertices(graph, mapped.communities)   
  graph.simplified <- simplify(graph.contracted,remove.multiple = FALSE ,remove.loops = FALSE)
  graph <- graph.simplified
  
  #Redefining attributes
  community.df <- data.frame(community = community, id = 1:n)
  #grouping all the members with same community
  community.grouped <- tapply(community.df$id, community.df$community, c)
  #Taking the sum of the attributes for the meta vertex
  mat <- sapply(community.grouped, function(x) colMeans(attr[x,]))    
  df2 <- as.data.frame(t(mat))
  row.names(df2) <- 1:nrow(df2)                                 
  attr <- df2
  
  #Updating the final answer
  ans2 <- list()
  names(community.grouped) <- 1:length(community.grouped)
  for(i in 1:length(community.grouped))
  {
    ans2[[i]] <- vector()
    for(j in 1:length(community.grouped[[i]]))
    {
      ans2[[i]] <- c(ans2[[i]], ans[[community.grouped[[i]][j]]])
    }
  }
  ans <- ans2
  filepath <- "~/Documents/06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection/communities_0.txt"
  if(compare(prev, community)$result == TRUE || length(ans)<6)
  {
    print("Phase 2 done!");
    flag = TRUE;
    for(x in ans)
    {
      if(flag)
      {
        write(x, file = filepath,sep = ",", ncolumns = 1000)    
        flag <- FALSE
      }
      else
      {
        write(x, file = filepath, append=TRUE, sep = ",", ncolumns = 1000)    
      }
    }
    break
  }
}

#Auxillary functions for distance computation
neu <- function(C,graph,x,community,mod)
{
  temp <- community
  temp[x] <- C
  return(modularity(graph, membership=temp)-mod)
}

phase1 <- function(community, alpha)
{
  epoch <- 15
  prev <- community
  for(k in 1:epoch)
  {
    for(i in 1:n)
    {
      x <- i                                     # x is the element to be moved
      original.mod <- modularity(graph, membership = community) # prev modularity before transfereing
      community.unique <- unique(community)                   # Set of unique communities
      
      neu.gains <- sapply(community.unique,neu,graph,x,community,original.mod)
      attr.gains <- sapply(community.unique,sim,graph,x,community,original.mod)
        #for debugging
#       if(length(neu.gains) < 32)
#       {
#         print(x)
#         print(length(community.unique))
#         print(neu.gains)
#         print(attr.gains/length(community.unique))
#         print("----------------------------------------")
#       }
      composite.gain <- alpha*neu.gains+(1-alpha)*attr.gains/length(community.unique)
      community.maxgain <- which.max(composite.gain)                    #community which has the best gain
      community[x] <- community.unique[community.maxgain]               #Assign x to that community
    }
    if(compare(community, prev)$result == TRUE)      #Check if community is changed or not
    {
      print("Converged Successfully!")
      break
    }
    prev <- community
  }
  fc<-community
  colors <- rainbow(max(fc))
  plot(simplify(graph),vertex.color=colors[fc], layout=layout.fruchterman.reingold, vertex.label= NA)
  return(community)
}

sim <- function(C,graph,x,community,mod)
{
  C.members <- which(community == C)
  Cx.members <- which(community == community[x])
  s_other <- 0
  s_this <- 0
  l1 <- 0
  l2 <- 0
  for(m in 1:length(C.members))
  {
    if(C.members[m]!=x)
    {
      s_other <- s_other + similarity[x,C.members[m]]
      
    }
    l1 <- l1 + length(ans[[C.members[m]]])
    
  }
  for(m in 1:length(Cx.members))
  {
    if(Cx.members[m]!=x)
    {
      s_this <- s_this + similarity[x,Cx.members[m]]
    }
    l2 <- l2 + length(ans[[Cx.members[m]]])
  }
  return(s_other/l1 - s_this/l2)
}
