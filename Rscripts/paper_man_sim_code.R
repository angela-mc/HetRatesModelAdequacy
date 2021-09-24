###################
## SIMULATION CODE#
###################

# (1) a single, internal branch shift not passed to descendants
# function takes a phylogeny, a rate-change magnitude, and how many branches change rates (no_events, here always 1)

scenario1 <- function (tree, magnitude,no_events)
{#get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es) #  one branch across all present analyses 
  demo_branch<- as.numeric(sample(int_br,no_events))
  
  # identify number of descendants for each of the int branches (default one internal branch)
  node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
  desc_no<-numeric() # number desc 
  desc<-list() # descendants
  for(i in 1:length(node_br)) # if one demo_branch, node_br will have one element, so the for goes from 1 to 1
  {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
   desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  
  #condition
  while(any(desc_no<4) || length(unique(unlist(desc)))!= length(unlist(desc)) ) #for 1+ shifts - do not have same descendants i.e. descendants do not repeat themselves
  {demo_branch<- as.numeric(sample(int_br,no_events))
   node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
   desc_no<-numeric()
   desc<-list()
   for(i in 1:length(node_br))
   {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
    desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  }
  
  # re-scale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[demo_branch] <- tree$edge.length[demo_branch]*magnitude
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[demo_branch]<- "red"
  edgecol[-demo_branch]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  return(list("tree"=tree, "demo_branch"=demo_branch))  
}

## apply as:
## scenario1 (tree, 10, 1) # rate change of x10, on 1 branch (no_events)


# (2) a clade event, in which all members of a particular group record a change in the rate of evolution 
# function takes a phylogeny, a rate-change magnitude, and how many clades change rates (no_events, here always 1)

scenario2 <- function (tree, magnitude,no_events) 
{#get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es), one clade across all present analyses
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of branch i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  # condition: clade must have between 15 and 30 members (=no_branches)
  while(any(desc_no2<15) || any(desc_no2>30) || length(unique(unlist(desc2)))!= length(unlist(desc2)) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   desc2<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  }
  
  # re-scale tree, by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  for(i in 1:length(node_br2)) # need for in case of multiple events
    tree$edge.length[as.numeric(names(desc2[[i]]))] <- tree$edge.length[as.numeric(names(desc2[[i]]))]*magnitude
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[as.numeric(names(desc2[[1]]))]<- "red"
  edgecol[-as.numeric(names(desc2[[1]]))]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  return(list("tree"=tree, "demo_branch2" = demo_branch2)) 
}

## apply as:
## scenario2 (tree, 10, 1) # rate-shifts on x10, on all members in 1 (no_events) clade



# (3) rate shifts on terminal branches
# function takes a phylogeny, a rate-change magnitude, and number of events (it's here just because it's in all rest of fxns, makes life easy when looping to all have same arguments)

scenario3 <- function (tree, magnitude,no_events) 
{ # get terminal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  
  #sample %no_branches & sort
  no_branches <- sample(15:30,1)
  modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
  
  #check condition: branches must not form a clade (i.e. consecutive tips)
  result <- rle(diff(modify_br)) #  values of differences & their frequency
  subset(as.list(result)$lengths, as.list(result)$values==1) -> condition # subset the ones with differences of 1 (i.e. the consecutive ones)
  # keep re-sampling until condition is OK
  while(all(condition>=3)) 
  {modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
   result <- rle(diff(modify_br))
   subset(as.list(result)$lengths, as.list(result)$values==1) -> condition
  }
  
  # rescale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[modify_br] <- tree$edge.length[modify_br]*magnitude
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[modify_br]<- "red"
  edgecol[-modify_br]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  
  return (list("tree"=tree, "modify_br"=modify_br)) 
}

## apply as: 
## scenario3 (tree, 10,1) # rate-shifts of x10 



# combinations
##############

# determines first X generations of node 'nodu' (here maxim = number of generations)
gener <- function(tree, nodu, start, maxim) #for now, it will also include the node about which you are asking, but it doesn't matter
{ rownames(tree$edge)<- 1:length(tree$edge.length)
  not_allowed<-getDescendants(tree,nodu)
  
  lngth=2
  if (length(not_allowed)<lngth)
    lngth=length(not_allowed)
  returned_value= nodu # a list of not_allowed nodes i.e. first descendants
  if (start < maxim) #get only the first maxim generations
    for(i in 1:lngth)
      if (not_allowed[i] != nodu) #so it won't put the same node many times
        returned_value = c(returned_value,gener(tree, not_allowed[i],start+1,maxim )) #get the children of the children 
  return (returned_value)
}

## apply as: 
## gener(tree, 15,0,4) # i.e first 4 generations of node=15
## as.numeric(names(gener(tree,15,0,4))) # gets the first X generations of branches for the target node


# (1+2) single internal burst + clade event
scenario12 <- function (tree, magnitude,no_events)
{ # do scenario 1 (single internal burst)
  #get int branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es) #  one branch across all present analyses
  demo_branch<- as.numeric(sample(int_br,no_events))
  
  # identify number of descendants for each of the int branches (default one branch)
  node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
  desc_no<-numeric() 
  desc<-list() # descendents
  for(i in 1:length(node_br)) # if 1 demo_branch, node_br will have one element, so the for goes from 1 to 1
  {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
   desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  
  #condition
  while(any(desc_no<4) || length(unique(unlist(desc)))!= length(unlist(desc)) ) #for 1+ shifts - do not have same descendants i.e. descendants do not repeat themselves
  {demo_branch<- as.numeric(sample(int_br,no_events))
   node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
   desc_no<-numeric()
   desc<-list()
   for(i in 1:length(node_br))
   {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
    desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  }
  
  # re-scale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[demo_branch] <- tree$edge.length[demo_branch]*magnitude
  
  
  # do 2, with the condition that:
  # you are not allowed to choose demo_branch2 from the ones corresponding to first x generations of (1)
  # simpler: demo_br 1 cannot be descendent of demo_br 2, because all descendants of demo_br 2 will be modified
  
  #get no_events random branch(es), default one clade event
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of demo_branch2 i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  #this only works for no_events =1, which is what I used across all analyses; for no_events >1, combination functions must be re-written
  # conditions 
  
  while(any(desc_no2<15) || any(desc_no2>30) || length(unique(unlist(desc2)))!= length(unlist(desc2)) || demo_branch2%in% as.numeric(names(gener(tree,node_br,0,4))) || demo_branch%in% names(desc2[[1]]) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   des2c<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  }
  
  # re-scale tree, by magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  for(i in 1:length(node_br2))
    tree$edge.length[as.numeric(names(desc2[[i]]))] <- tree$edge.length[as.numeric(names(desc2[[i]]))]*magnitude
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[as.numeric(names(desc2[[1]]))]<- "red"
  edgecol[demo_branch]<- "green"
  edgecol[-c(as.numeric(names(desc2[[1]])), demo_branch)]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  
  return(list("tree"=tree, "demo_branch"= demo_branch, "demo_branch2" = demo_branch2))
}

## apply as:
## scenario12(tree,10,1) # 1 clade and 1 sg-branch have x10 rate-shifts


# (1+3) single internal burst + rate changes on terminal branches
# do (1) first, then (3) with the condition that branches are not from the ones x descendants of demo_branch

scenario13 <- function (tree, magnitude,no_events)
{ # do scenario 1 - single internal burst
  #get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es) #  one branch across all present analyses 
  demo_branch<- as.numeric(sample(int_br,no_events))
  
  # identify number of descendants for each of the int branches (default one branch)
  node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
  desc_no<-numeric() 
  desc<-list() # descendants
  for(i in 1:length(node_br)) # if 1 demo_branch, node_br will have one element, so the for goes from 1 to 1
  {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
   desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  
  #condition
  while(any(desc_no<4) || length(unique(unlist(desc)))!= length(unlist(desc)) ) #for 1+ shifts - do not have same descendants i.e. descendants do not repeat themselves
  {demo_branch<- as.numeric(sample(int_br,no_events))
   node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
   desc_no<-numeric()
   desc<-list()
   for(i in 1:length(node_br))
   {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
    desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  }
  
  # re-scale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[demo_branch] <- tree$edge.length[demo_branch]*magnitude
  
  # # do shifts on terminal branches, with condition
  
  #sample %no_branches & sort
  no_branches <- sample(15:30,1)
  modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
  
  #check condition - tips are not consecutive
  result <- rle(diff(modify_br)) #  values of differences & their frequency
  subset(as.list(result)$lengths, as.list(result)$values==1) -> condition # subset the ones with differences of 1 (i.e. the consecutive tips)
  
  
  # keep re-sampling until condition is OK & no overlap between scenarios
  while(all(condition>=3) || any(modify_br%in% as.numeric(names(gener(tree,node_br,0,4))))  ) 
  {modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
   result <- rle(diff(modify_br))
   subset(as.list(result)$lengths, as.list(result)$values==1) -> condition
  }
  
  # rescale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[modify_br] <- tree$edge.length[modify_br]*magnitude
  
  # plot colour coded
  
  edgecol<-tree$edge.length
  edgecol[modify_br]<- "red"
  edgecol[demo_branch]<- "green"
  edgecol[-c(modify_br, demo_branch)]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  return(list("tree"=tree, "demo_branch"=demo_branch, "modify_br"=modify_br))
}

## apply as:
## scenario13(tree,10,1) # terminal branches and a single branch have shifts of x10


## (2+3) clade event +  changes on terminal branches 
#do (2) first, and then realistically 3 should not be in clade 2 (or in demo branch and descendants for simplicity)

scenario23<- function(tree, magnitude,no_events)
{# do 2
  #get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es), one clade event across all present analyses
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of branch i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  # condition: clade has 15-30 members (=no_branches)
  while(any(desc_no2<15) || any(desc_no2>30)|| length(unique(unlist(desc2)))!= length(unlist(desc2)) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   desc2<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  }
  
  # re-scale tree, by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  for(i in 1:length(node_br2)) # needed for in case of multiple events - here not the case
    tree$edge.length[as.numeric(names(desc2[[i]]))] <- tree$edge.length[as.numeric(names(desc2[[i]]))]*magnitude
  
  
  # do 3, with condition
  no_branches<-sample(15:30,1)
  modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
  
  #check condition: tips are not consecutive
  result <- rle(diff(modify_br)) #  values of differences & their frequency
  subset(as.list(result)$lengths, as.list(result)$values==1) -> condition # subset the ones with differences of 1 (i.e. the consecutive ones)
  # keep re-sampling until condition is OK
  
  #this only works for no_events =1 (what I used across all analyses); if more events, one should modify condition
  while(all(condition>=3) || any(modify_br%in% names(desc2[[1]]))  ) 
  {modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
   result <- rle(diff(modify_br))
   subset(as.list(result)$lengths, as.list(result)$values==1) -> condition
  }
  
  # rescale tree, by magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[modify_br] <- tree$edge.length[modify_br]*magnitude
  
  # # plot colour coded
  
  edgecol<-tree$edge.length
  edgecol[modify_br]<- "red"
  edgecol[as.numeric(names(desc2[[1]]))]<- "green"
  edgecol[-c(modify_br, as.numeric(names(desc2[[1]])))]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  
  return(list("tree"=tree, "demo_branch2"=demo_branch2, "modify_br"=modify_br))
}

## apply as:
## scenario23(tree,10,1) # one clade and terminal branches have x10 rate shifts


# (1+2+3) single-lineage burst + clade event + changes on terminal branches
# (1)
# then (2), condition not within the first 3 generations & ancestors
# then (3), condition outside (2) & not within 3 generations from (1)

scenario123 <- function(tree, magnitude,no_events)
{# do scenario 1 - single internal branch shift
  #get int branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es) #  one branch across all present analyses
  demo_branch<- as.numeric(sample(int_br,no_events))
  
  # identify number of descendants for each of the int branches (default one branch)
  node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
  desc_no<-numeric() 
  desc<-list() # descendants
  for(i in 1:length(node_br)) # if 1 demo_branch, node_br will have one element, so the for goes from 1 to 1
  {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
   desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  
  #condition
  while(any(desc_no<4) || length(unique(unlist(desc)))!= length(unlist(desc)) ) #for 1+ shifts - do not have same descendants i.e. descendants do not repeat themselves
  {demo_branch<- as.numeric(sample(int_br,no_events))
   node_br<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch)[,2] # node at the end of demo_branch
   desc_no<-numeric()
   desc<-list()
   for(i in 1:length(node_br))
   {desc_no[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2])))
    desc[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch[i])[,2]))}
  }
  
  # re-scale tree by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[demo_branch] <- tree$edge.length[demo_branch]*magnitude
  
  
  # do scenario 2 (clade event), with the condition that:
  # you are not allowed to choose demo_branch2 from the ones corresponding to first x generations of (1)
  # simpler: demo_br (1) cannot be descendent of demo_br (2), because all descendants of demo_br (2) will be modified
  
  #get no_events random branch(es), default one clade event
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of demo_branch2 i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  #this only works for no_events =1, which is what I used across all analyses; for more events, one should modify code
  # conditions 
  # clade has 15-30 memebers(no branches)
  
  while(any(desc_no2<15) || any(desc_no2>30) || length(unique(unlist(desc2)))!= length(unlist(desc2)) || demo_branch2%in% as.numeric(names(gener(tree,node_br,0,4))) || demo_branch%in% names(desc2[[1]]) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   des2c<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  }
  
  # re-scale tree, by magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  for(i in 1:length(node_br2))
    tree$edge.length[as.numeric(names(desc2[[i]]))] <- tree$edge.length[as.numeric(names(desc2[[i]]))]*magnitude
  
  
  # do scenario 3 (rate changes at tips), with conditions - not in first X gen from (1) & not in the members of (2)
  no_branches<-sample(15:30,1)
  modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
  
  #check condition: not consecutive tips
  result <- rle(diff(modify_br)) #  values of differences & their frequency
  subset(as.list(result)$lengths, as.list(result)$values==1) -> condition # subset the ones with differences of 1 (i.e. the consecutive ones)
  # keep re-sampling until condition is OK
  
  #this only works for no_events =1; if more, one should modify condition
  while(all(condition>=3) || any(modify_br%in% as.numeric(names(gener(tree,node_br,0,4)))) || any(modify_br%in% names(desc2[[1]]))) 
  {modify_br<- sort(sample(t_br, round((no_branches/100)*length(tree$tip.label))))
   result <- rle(diff(modify_br))
   subset(as.list(result)$lengths, as.list(result)$values==1) -> condition
  }
  
  # rescale tree, by magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  tree$edge.length[modify_br] <- tree$edge.length[modify_br]*magnitude
  
  # # plot colour coded
  
  edgecol<-tree$edge.length
  edgecol[demo_branch]<- "blue"
  edgecol[modify_br]<- "red"
  edgecol[as.numeric(names(desc2[[1]]))]<- "green"
  edgecol[-c(demo_branch, modify_br, as.numeric(names(desc2[[1]])))]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  
  return(list("tree"=tree, "demo_branch"=demo_branch, "demo_branch2"=demo_branch2, "modify_br"=modify_br))
}

## apply as: 
## scenario123(tree, 10,1) # an internal single branch, one clade and terminal branches have x10 rate-shifts 


## (4) a constant rate-deceleration process from the root to tips 
build_root_process<-function(rate_decr, size, number_sim)
{# simulate original trees
  sim_trees<- sim.bd.taxa(n=size,numbsim=number_sim,lambda=1,mu=0,complete=FALSE)
  names(sim_trees)<-list()
  rescaled_trees <- list()
  names(rescaled_trees)<-list()
  sim_data<-list()
  names(sim_data)<-list()
  
  for(i in 1:number_sim)
  {#standardize depth = 1
    sim_trees[[i]]$edge.length<- sim_trees[[i]]$edge.length * ( 1/ (  diag(vcv.phylo(sim_trees[[i]]))[1] ) )
    # create nex file
    write.nexus (sim_trees[[i]], file=paste("root_process","_", rate_decr,"_",i,"nexus.nex", sep=""))
    names(sim_trees)[[i]] <- paste("root_process",rate_decr,i, sep="_")
    
    # rescale tree
    rescaled_trees[[i]]<-rescale(sim_trees[[i]], model=c("EB"), a=log(rate_decr))
    names(rescaled_trees)[[i]] <- paste("root_process", rate_decr,i, sep="_")
    
    # simulate data
    sim_data[[i]]<- sim.char(rescaled_trees[[i]], 1, 1, root=100)
    names(sim_data)[[i]]<- paste("root_process", rate_decr,i, sep="_")
    #create txt file
    write.table(sim_data[[i]], file=paste("root_process", "_", rate_decr,"_",i,".txt", sep=""), quote=FALSE, col.names = FALSE, sep="\t")  
  }
  save(sim_trees, file=paste("original","_","root_process","_", rate_decr, ".rda", sep=""))
  save(rescaled_trees, file=paste("rescaled","_","root_process","_", rate_decr, ".rda", sep=""))
  save(sim_data, file=paste("data","_","root_process", "_",rate_decr, ".rda", sep="") )
  names(sim_trees)-> list_names
  save(list_names, file=paste("names", "_", "root_process", "_", rate_decr, ".rda", sep=""))
}


## (5) a case when single clade goes through an initial increase in the rate of evolution (x5) followed by a constant rate-decay

## (1) extract clade
my.extract.clade<-function(phy,node,root.edge=0){
  #requires ape and phylobase be loaded
  library(phylobase)
  phy.phylo4<-phylo4(phy)
  save<-descendants(phy.phylo4,node)
  extracted<-drop.tip(phy,(1:length(phy$tip.label))[(1:length(phy$tip.label))%in%save==F],root.edge=root.edge)
  return(extracted)
}

## (2) simulate scenario
clade_process <- function (tree, magnitude,no_events) 
{#get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es), one clade event through-out the analyses
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of branch i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  # condition: 15-30 members (=no_branches)
  while(any(desc_no2<15) || any(desc_no2>30) || length(unique(unlist(desc2)))!= length(unlist(desc2)) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   desc2<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  }
  
  # re-scale clade, by magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  new_tree <- my.extract.clade(tree,node_br2)
  
  new_tree$edge.length<-new_tree$edge.length*5
  rescaled_new<-rescale(new_tree, model = c("EB"), a=log(magnitude))
  scaled<-rescaled_new$edge.length/new_tree$edge.length
  for(i in 1:length(node_br2)) # need for in case of multiple events; not the case in the present analyses though
    tree$edge.length[ sort(as.numeric(names(desc2[[i]]))) ] <-rescaled_new$edge.length
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[as.numeric(names(desc2[[1]]))]<- "red"
  edgecol[-as.numeric(names(desc2[[1]]))]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
  return(list("tree"=tree, "demo_branch2" = demo_branch2)) 
}


## build trees & associated trait-data for scenarios 1,2,3 & combinations
## function takes: scenario_function, shift magnitude, size of tree, number of simulations and name of scenario (character)

build_trees<-function(scenario, magnitude, size, number_sim, scenario_name)
{# simulate input trees
  sim_trees<- sim.bd.taxa(n=size,numbsim=number_sim,lambda=1,mu=0,complete=FALSE)
  names(sim_trees)<-list()
  rescaled_trees <- list()
  names(rescaled_trees)<-list()
  sim_data<-list()
  names(sim_data)<-list()
  
  for(i in 1:number_sim)
  {#standardize depth = 1
    sim_trees[[i]]$edge.length<- sim_trees[[i]]$edge.length * ( 1/ (  diag(vcv.phylo(sim_trees[[i]]))[1] ) )
    # create nex file
    write.nexus (sim_trees[[i]], file=paste(scenario_name,"_", magnitude,"_",i,"nexus.nex", sep=""))
    names(sim_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
    
    # rescale tree
    rescaled_trees[[i]]<-scenario(sim_trees[[i]], magnitude,1) 
    #rescaled_trees[[i]][[1]] = trees, rest is branches modified
    names(rescaled_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
    
    # simulate data
    sim_data[[i]]<- sim.char(rescaled_trees[[i]][[1]], 1, 1, root=100)
    names(sim_data)[[i]]<- paste(scenario_name, magnitude,i, sep="_")
    #create txt file
    write.table(sim_data[[i]], file=paste(scenario_name, "_", magnitude,"_",i,".txt", sep=""), quote=FALSE, col.names = FALSE, sep="\t")  
  }
  save(sim_trees, file=paste("original","_",scenario_name,"_", magnitude, ".rda", sep=""))
  save(rescaled_trees, file=paste("rescaled","_",scenario_name,"_", magnitude, ".rda", sep=""))
  save(sim_data, file=paste("data","_",scenario_name, "_",magnitude, ".rda", sep="") )
  names(sim_trees)-> list_names
  save(list_names, file=paste("names", "_", scenario_name, "_", magnitude, ".rda", sep=""))
}

## apply as: 
## build_trees(scenario1, 2, size=100, number_sim=100, "scenario1") # create 100 trees, size=100 tips, heterogeneity scenario: single, internal branch shift, rate-shift = 2

##---------------------------------------------------------------------------


###  same functions work for trees of sizes: 25,50,100,200; however, in these analyses the number of branches that change the rate was set to 10-15
### problem with a clade event is that: in a 25-tips tree, might not find a clade that has 10-15 members, therefore it needs a changed function

## clade-event (scenario2) for the analyses involving differed sized trees
## additional comment - might not find a clade 10-15branches in a tree of 25 species
## condition that if it doesn't find one --> re-simulate tree

scenario2_25 <- function (tree, magnitude,no_events) 
{#get internal branches
  rownames(tree$edge)<- 1:length(tree$edge.length)
  t_br<-as.numeric(rownames(subset(tree$edge, tree$edge[,2]<=length(tree$tip.label))))
  int_br<-(1:length(tree$edge.length))[-t_br]
  
  #get no_events random branch(es), default should be one clade event
  demo_branch2<- as.numeric(sample(int_br,no_events))
  
  # descendants of branch i.e. clade to be modified
  node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
  desc_no2<-numeric()
  desc2<-list()
  for(i in 1:length(node_br2))
  {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
   desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
  
  # condition: 10-15 members in a clade (=no_branches)
  
  all_branches<-list()
  all_branches[[1]]<-demo_branch2
 
  while( (any(desc_no2<10) || any(desc_no2>15) || length(unique(unlist(desc2)))!= length(unlist(desc2))) && (length(all_branches)<length(tree$edge.length)) )
  {demo_branch2<- as.numeric(sample(int_br,no_events))
   node_br2<-subset(tree$edge, names(tree$edge[,2])%in%demo_branch2)[,2]
   desc_no2<-numeric()
   desc2<-list()
   for(i in 1:length(node_br2))
   {desc_no2[i]<-length(getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2])))
    desc2[[i]]<- getDescendants(tree, (subset(tree$edge, names(tree$edge[,2])%in%demo_branch2[i])[,2]))}
   
   all_branches<-c(all_branches, demo_branch2)
  }
  
 if(length(all_branches)<length(tree$edge.length)) # if there is a good demo_branch2
 {# re-scale tree, by *magnitude
  names(tree$edge.length) <- 1:length(tree$edge.length)
  for(i in 1:length(node_br2)) # need for in case of multiple events
    tree$edge.length[as.numeric(names(desc2[[i]]))] <- tree$edge.length[as.numeric(names(desc2[[i]]))]*magnitude
  
  # colour coded plot
  edgecol<-tree$edge.length
  edgecol[as.numeric(names(desc2[[1]]))]<- "red"
  edgecol[-as.numeric(names(desc2[[1]]))]<- "black"
  
  plot(tree,cex=0.5, edge.color = edgecol)
}
  if(length(all_branches)==length(tree$edge.length)) demo_branch2<-NULL
   
  return(list("tree"=tree, "demo_branch2" = demo_branch2)) 
}


### build trees --> incorporates the problem with clade size and scenario 2 (clade event) for 25 species
build_trees<-function(scenario, magnitude, size, number_sim, scenario_name)
{# simulate original trees
  sim_trees<- sim.bd.taxa(n=size,numbsim=number_sim,lambda=1,mu=0,complete=FALSE)
  names(sim_trees)<-list()
  rescaled_trees <- list()
  names(rescaled_trees)<-list()
  sim_data<-list()
  names(sim_data)<-list()
  
  for(i in 1:number_sim)
  {#standardize depth = 1
    sim_trees[[i]]$edge.length<- sim_trees[[i]]$edge.length * ( 1/ (  diag(vcv.phylo(sim_trees[[i]]))[1] ) )
    # create nex file
    write.nexus (sim_trees[[i]], file=paste(scenario_name,"_", magnitude,"_",i,"nexus.nex", sep=""))
    names(sim_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
    
    if(scenario_name%in%c("scenario2") && size==25) # for clade-event, in trees of 25 species
      { # rescale tree
        rescaled_trees[[i]]<-scenario2_25(sim_trees[[i]], magnitude,1) 
        #rescaled_trees[[i]][[1]] = rescaled trees
        names(rescaled_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
      while(is.null(rescaled_trees[[i]]$demo_branch))
        { sim_trees[[i]]<- sim.bd.taxa(n=size,numbsim=1,lambda=1,mu=0,complete=FALSE)[[1]] # replace tree with another simulated tree
          sim_trees[[i]]$edge.length<- sim_trees[[i]]$edge.length * ( 1/ (  diag(vcv.phylo(sim_trees[[i]]))[1] ) )
          # re-write tree: create nex file 
          write.nexus (sim_trees[[i]], file=paste(scenario_name,"_", magnitude,"_",i,"nexus.nex", sep=""))
          names(sim_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
        
          rescaled_trees[[i]]<-scenario2_25(sim_trees[[i]], magnitude,1) 
          #rescaled_trees[[i]][[1]] = rescaled 
          names(rescaled_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")
          #print(paste(scenario_name, magnitude,i, sep="_"))
      }}
    
    ## if not - proceed as normal
    # rescale tree
    rescaled_trees[[i]]<-scenario(sim_trees[[i]], magnitude,1) 
    #rescaled_trees[[i]][[1]] = trees, rest is branches modified
    names(rescaled_trees)[[i]] <- paste(scenario_name, magnitude,i, sep="_")    
    
    # simulate data
    sim_data[[i]]<- sim.char(rescaled_trees[[i]][[1]], 1, 1, root=100)
    names(sim_data)[[i]]<- paste(scenario_name, magnitude,i, sep="_")
    #create txt file
    write.table(sim_data[[i]], file=paste(scenario_name, "_", magnitude,"_",i,".txt", sep=""), quote=FALSE, col.names = FALSE, sep="\t")  
  }
  save(sim_trees, file=paste("original","_",scenario_name,"_", magnitude, ".rda", sep=""))
  save(rescaled_trees, file=paste("rescaled","_",scenario_name,"_", magnitude, ".rda", sep=""))
  save(sim_data, file=paste("data","_",scenario_name, "_",magnitude, ".rda", sep="") )
  names(sim_trees)-> list_names
  save(list_names, file=paste("names", "_", scenario_name, "_", magnitude, ".rda", sep=""))
}
