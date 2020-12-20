## -----------------------------------------------------------------------------
library(igraph)

## -----------------------------------------------------------------------------
# Prim function
Prim <- function(vn, adjacency_matrix, myl){
  link_p <- function(x,y,w){#构建邻接矩阵与子节点矩阵
  adjacency_matrix_p[x,y]<<-w
  adjacency_matrix_p[y,x]<<-w
  }
  adjacency_matrix_p=matrix(0,vn,vn,byrow=T,dimnames=list(letters[1:vn],letters[1:vn])) 
  a<-1
  visited <- c(0*(1:vn))
  visited[a]=1
  myg <- graph.adjacency(adjacency_matrix,mode = "undirected",weighted = TRUE)
  for(o in 1:vn-1 ){
    for(i in order(E(myg)$weight)){
      i=as.numeric(i)
      v1<-as.numeric(charToRaw(get.edgelist(myg)[i,1]))-as.numeric(charToRaw('a'))+1
      v2<-as.numeric(charToRaw(get.edgelist(myg)[i,2]))-as.numeric(charToRaw('a'))+1
      if(visited[v1]==1&&visited[v2]==0){
        link_p(v1,v2,E(myg)[i]$weight)
        visited[v2]=1
        E(myg)[i]$color="red"
        plot.igraph(myg,layout=myl,main="Prim",edge.label = E(myg)$weight)
        Sys.sleep(1)
        break
      }
      if(visited[v2]==1&&visited[v1]==0){
        link_p(v1,v2,E(myg)[i]$weight)
        visited[v1]=1
        E(myg)[i]$color="red"
        plot.igraph(myg,layout=myl,main="Prim",edge.label = E(myg)$weight)
        Sys.sleep(1)
        break
      }
    }
  }}

## -----------------------------------------------------------------------------
vn <- 12
a <- c(0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0,
 2, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0,
 0, 3, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0,
 3, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 0,
 0, 1, 0, 0, 4, 0, 3, 0, 0, 2, 0, 0,
 0, 0, 2, 0, 0, 3, 0, 3, 0, 0, 4, 3,
 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 4, 0, 0, 0, 0, 3, 0, 0,
 0, 0, 0, 0, 0, 2, 0, 0, 3, 0, 3, 0,
 0, 0, 0, 0, 0, 0, 4, 0, 0, 3, 0, 1,
 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0)
adjacency_matrix <- matrix(data = a, nrow = 12, ncol = 12, dimnames=list(letters[1:vn],letters[1:vn]))
b <- c(0,1,2,3,0,1,2,3,0,1,2,3,2,2,2,2,1,1,1,1,0,0,0,0)
myl <- matrix(data = b, nrow = 12, ncol = 2)
rm(a)


Prim(vn, adjacency_matrix, myl)

## -----------------------------------------------------------------------------
Kruskal <- function(vn, adjacency_matrix, myl){
  bindset <- function(x,y){
  same_set[x,y]<<-1
  same_set[y,x]<<-1
  for(i in 1:vn){
  i=as.numeric(i)
  if(same_set[y,i]==1&&i!=x&&i!=y&&same_set[x,i]==0) bindset(i,x)
  if(same_set[x,i]==1&&i!=y&&i!=x&&same_set[y,i]==0) bindset(i,y)
  }
}
  myg <- graph.adjacency(adjacency_matrix,mode = "undirected",weighted = TRUE)
  cnt<-0
  link_k <- function(x,y,w){#构建邻接矩阵与子节点矩阵
  adjacency_matrix_k[x,y]<<-w
  adjacency_matrix_k[y,x]<<-w
  }
  adjacency_matrix_k=matrix(0,vn,vn,byrow=T,dimnames=list(letters[1:vn],letters[1:vn])) #创建k子树邻接矩阵
  for(i in order(E(myg)$weight)){
    i=as.numeric(i)
    v1<-as.numeric(charToRaw(get.edgelist(myg)[i,1]))-as.numeric(charToRaw('a'))+1
    v2<-as.numeric(charToRaw(get.edgelist(myg)[i,2]))-as.numeric(charToRaw('a'))+1
    if(myset[v1]==0&&myset[v2]==0){
      setcnt<<-setcnt+1
      myset[v1]<<-setcnt
      myset[v2]<<-setcnt
      link_k(v1,v2,E(myg)[i]$weight)
      E(myg)[i]$color="red"
      plot.igraph(myg,layout=myl,main="Kruskal",edge.label = E(myg)$weight)
      #Sys.sleep(1)
      cnt=cnt+1
      if(cnt==vn-1) break
    }
    else{ if(myset[v1]==0||myset[v2]==0){
      if(myset[v2]==0) myset[v2]=myset[v1]
      else myset[v1]=myset[v2]
      link_k(v1,v2,E(myg)[i]$weight)
      E(myg)[i]$color="red"
      plot.igraph(myg,layout=myl,main="Kruskal",edge.label = E(myg)$weight)
      Sys.sleep(1)
      cnt=cnt+1
      if(cnt==vn-1) break
    }
    else if(same_set[myset[v1],myset[v2]]==0){
      bindset(myset[v1],myset[v2])
      link_k(v1,v2,E(myg)[i]$weight)
      E(myg)[i]$color="red"
      plot.igraph(myg,layout=myl,main="Kruskal",edge.label = E(myg)$weight)
      Sys.sleep(1)
      cnt=cnt+1
      if(cnt==vn-1) break
    }
    }
  }
}

## -----------------------------------------------------------------------------
vn <- 12
a <- c(0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0,
 2, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0,
 0, 3, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0,
 3, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 0,
 0, 1, 0, 0, 4, 0, 3, 0, 0, 2, 0, 0,
 0, 0, 2, 0, 0, 3, 0, 3, 0, 0, 4, 3,
 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 4, 0, 0, 0, 0, 3, 0, 0,
 0, 0, 0, 0, 0, 2, 0, 0, 3, 0, 3, 0,
 0, 0, 0, 0, 0, 0, 4, 0, 0, 3, 0, 1,
 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0)
adjacency_matrix <- matrix(data = a, nrow = 12, ncol = 12, dimnames=list(letters[1:vn],letters[1:vn]))
b <- c(0,1,2,3,0,1,2,3,0,1,2,3,2,2,2,2,1,1,1,1,0,0,0,0)
myl <- matrix(data = b, nrow = 12, ncol = 2)
rm(a)
myset<-c(0*(1:vn))
same_set=matrix(0,vn,vn,byrow=T)
setcnt<-0
for(i in 1:vn) same_set[i,i]=1

Kruskal(vn, adjacency_matrix, myl)

