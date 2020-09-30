# I commented out all the commands to output pdf and csv files
# This code should output tables with statistics for minifish repertoire similarity networks as
# well as some barplots and simulation of node removal

library(igraph)
library(poweRlaw)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# Load data, collapse into aminoacid CDR3 sequences
  load(file = "MF.list")
  MF.x.list  = MF.list
  MF.fish.names = c("Fish 2","Fish 3","Fish 4","Fish 5")
  MF.gene.names = c("tcra","tcrb","tcrg","tcrd","igh")
  names(MF.x.list) = MF.gene.names
  for (i in 1:length(MF.x.list))names(MF.x.list[[i]]) = MF.fish.names
# Collapse allelic variants of Js and make a dataframe
  MF.Ja.alleles.names = c(1,6,10,12,21,30,31,34,39,42,43,45,47,48,50,51,55,60) 
  for(i in 1:length(MF.x.list$tcra)){
    ja.alleles = MF.Ja.alleles.names[match(MF.x.list$tcra[[i]]$J.gene,62:79)]
    MF.x.list$tcra[[i]]$J.gene[!is.na(ja.alleles)] = ja.alleles[!is.na(ja.alleles)]
  }
  MF.x.df = lapply(MF.x.list,function(x)do.call(rbind,x))
  MF.x.df = lapply(MF.x.df,function(x)cbind(x,fish=sub("\\..*","",rownames(x))))

  MF.aa.f = lapply(1:4,function(f)lapply(1:5,function(i)data.frame(MF.x.df[[i]]%>%filter(fish == MF.fish.names[[f]], L%%3==0, !grepl("\\*",CDR3.amino.acid.sequence))%>%group_by(CDR3.amino.acid.sequence)%>%summarize(UMI = sum(Umi.count),clones = n(),vs.u=n_distinct(V.gene),js.u=n_distinct(J.gene),pub = max(pub)))))
  names(MF.aa.f) = MF.fish.names
  
# make distance matrix and graph g.Ln.f[[fish]][[connect by distance d]][[gene]]
  sD.MF.fish.f = lapply(1:4,function(f)lapply(1:5,function(j)stringDist(MF.aa.f[[f]][[j]]$CDR3.amino.acid.sequence)))
  mat.MF.f =lapply(1:4,function(f)lapply(sD.MF.fish.f[[f]],as.matrix))
  g.Ln.f = lapply(1:4,function(f)lapply(1:12,function(n)lapply(mat.MF.f[[f]],function(x)graph_from_adjacency_matrix(x==n,mode = "undirected"))))
  
# add information about V, J and graph statistics
  for(f in 1:4){
    for (i in 1:5){
      MF.aa.f[[f]][[i]]$V = 0
      MF.aa.f[[f]][[i]]$V[MF.aa.f[[f]][[i]]$vs.u ==1] = MF.x.list[[i]][[f]]$V.gene[match(MF.aa.f[[f]][[i]][MF.aa.f[[f]][[i]]$vs.u ==1,]$CDR3.amino.acid.sequence,MF.x.list[[i]][[f]]$CDR3.amino.acid.sequence)]
      MF.aa.f[[f]][[i]]$J = 0
      MF.aa.f[[f]][[i]]$J[MF.aa.f[[f]][[i]]$js.u ==1] = MF.x.list[[i]][[f]]$J.gene[match(MF.aa.f[[f]][[i]][MF.aa.f[[f]][[i]]$js.u ==1,]$CDR3.amino.acid.sequence,MF.x.list[[i]][[f]]$CDR3.amino.acid.sequence)]
      MF.aa.f[[f]][[i]]$membership = components(g.Ln.f[[f]][[1]][[i]])$membership
      MF.aa.f[[f]][[i]]$csize = components(g.Ln.f[[f]][[1]][[i]])$csize[components(g.Ln.f[[f]][[1]][[i]])$membership]
      MF.aa.f[[f]][[i]]$authority = authority_score(g.Ln.f[[f]][[1]][[i]])$vector
      MF.aa.f[[f]][[i]]$degree = degree(g.Ln.f[[f]][[1]][[i]])
      MF.aa.f[[f]][[i]]$betweenness = betweenness(g.Ln.f[[f]][[1]][[i]])
      MF.aa.f[[f]][[i]]$coreness = coreness(g.Ln.f[[f]][[1]][[i]])
    }
  }

#############
### PLOTS ###
#############

# Dataframe for plots  
  MF.aa.df.f = list()
  for(f in 1:4)names(MF.aa.f[[f]]) = factor(c('tcra','tcrb','tcrg','tcrd','igh'))
  for(f in 1:4)MF.aa.df.f[[f]]  = do.call(rbind,MF.aa.f[[f]])
  for(f in 1:4)MF.aa.df.f[[f]]$gen = factor(gsub("\\..*","",rownames(MF.aa.df.f[[f]])),levels = c('tcra','tcrb','tcrg','tcrd','igh'))
  MF.aa.df.f = do.call(rbind,lapply(1:4,function(f)cbind(MF.aa.df.f[[f]],fish = names(MF.aa.f)[f])))

# Number of distinct Js
  mean(data.frame(MF.aa.df.f%>%filter(csize>4)%>%group_by(gen,membership,fish)%>%summarise(n = n_distinct(J)))$n)
  mean(data.frame(MF.aa.df.f%>%filter(csize>4)%>%group_by(gen,membership,fish)%>%summarise(n = n_distinct(V)))$n)
  df.js = MF.aa.df.f%>%filter(csize>4)%>%group_by(gen,membership,fish)%>%summarise(n = n_distinct(J)) %>% group_by(fish,gen) %>%summarise(n_mean = mean(n))
  df.vs = MF.aa.df.f%>%filter(csize>4)%>%group_by(gen,membership,fish)%>%summarise(n = n_distinct(V)) %>% group_by(fish,gen) %>%summarise(n_mean = mean(n))

# Two tables of V usage and J usage per cluster
  data.frame(df.js %>% spread(fish,n_mean))
  data.frame(df.vs %>% spread(fish,n_mean))
  
  #  pdf(file = "Minifish paper//js.vs_per_cluster_plots.pdf",width = 4,height = 3)
  ggplot(df.js)   + geom_boxplot(aes(y=n_mean, x = gen,colour = gen)) + geom_point(aes(y=n_mean, x = gen)) + scale_colour_manual(values = c(gg_color_hue(4)[c(1,3,2,4)],"brown"))
  ggplot(df.vs)   + geom_boxplot(aes(y=n_mean, x = gen,colour = gen)) + geom_point(aes(y=n_mean, x = gen)) + scale_colour_manual(values = c(gg_color_hue(4)[c(1,3,2,4)],"brown"))
  #  dev.off()
  
  # Make csv files
  #write.csv2(data.frame(df.js %>% spread(fish,n_mean)),file = "Minifish paper/js_per_cluster.csv")
  #write.csv2(data.frame(df.vs %>% spread(fish,n_mean)),file = "Minifish paper/vs_per_cluster.csv")

# Cluster size barplot
  #pdf(file = "Minifish paper//cluster size barplot.pdf",width = 20,height = 15)
  ggplot(MF.aa.df.f%>% filter(fish == "Fish 5")) + geom_bar(mapping = aes(x =factor(csize)), fill = "steel blue") +  guides(fill=guide_legend(title="Publicity")) + theme(axis.text.x = element_text(angle = 90,size = 14,margin = margin(t = c(20))),axis.text.y = element_text(size = 14),axis.title = element_text(size = 20)) + xlab ("Cluster size") + ylab ("CDR3s per cluster size") + facet_wrap(~gen,scales = "free") + theme(strip.text = element_text(face = "italic",size = 15))
  #dev.off()


# Jitter plot
  
  #pdf(file = "Minifish/network jitter.pdf")
  ggplot(MF.aa.df.f %>% filter(fish == "Fish 5")) + geom_jitter(mapping = (aes(color= factor(pub),x =log10(csize),y = pub)),size = 0.01,width = 0.05) + guides(color=guide_legend(title="Publicity")) + scale_colour_manual(values = c("black","dark green", "blue","red")) + ylab("Publicity") + xlab ("log10 cluster size") + facet_wrap(~gen,nrow = 5) + theme(strip.text = element_text(face = "italic"))
  #dev.off()
  # source data
  jitter.plot.source = sapply(1:5,function(i)table(components(g.Ln.f[[4]][[1]][[i]])$csize[components(g.Ln.f[[4]][[1]][[i]])$membership],MF.aa.f[[4]][[i]]$pub))
  jitter.plot.source = do.call(rbind,lapply(1:5,function(i)rbind(c(MF.gene.names[[i]],"","",""),jitter.plot.source[[i]])))
  jitter.plot.source
  #write.csv2(jitter.plot.source,file = "Minifish paper/jitter.plot.source.csv")
  
# Publicity - Degree boxplot
  
  #pdf(file = "Minifish paper/degree boxplot.pdf")
  df.deg.pub = Reduce('+',lapply(MF.aa.f,function(x)sapply(x,function(y)tapply(y$degree,y$pub,mean))) )/4
  df.deg.pub
  #write.csv2(df.deg.pub,file = "Minifish paper/boxplot.deg.pub.source.csv")
  
  ggplot(MF.aa.df.f %>% filter(fish =="Fish 5")) + geom_boxplot(mapping = (aes(group= pub,y = degree,color = factor(pub)))) + 
    guides(color=guide_legend(title="Publicity")) + 
    scale_colour_manual(values = c("black","dark green", "blue","red")) +
    facet_wrap(~gen)  + theme(strip.text = element_text(face = "italic"))
  #dev.off()

# Size largest component
  lc=sapply(1:4,function(f)sapply(1:5,function(gen)max(sort(unique(MF.aa.f[[f]][[gen]]$csize))/length(MF.aa.f[[f]][[gen]]$csize))))
  sapply(1:4,function(f)sapply(1:5,function(gen)max(sort(unique(MF.aa.f[[f]][[gen]]$csize))/length(MF.aa.f[[f]][[gen]]$csize))))
  rownames(lc) = MF.gene.names
  colnames(lc) = MF.fish.names
  lc
  #write.csv2(lc,file = "Minifish paper/size of largest component.csv")
  
  
# Node removal

# Removing IgH for all genes
{

  p = list()
  deg.maxs = list()
  for (fish in 1:4){
    for (gen in 5:5){
      thr.pub = 2
      sims = 40
      remove = floor(sum(MF.aa.f[[fish]][[gen]]$pub>=thr.pub)*0.66)
      #remove = 20
      #remove = 200  
      g = g.Ln.f[[fish]][[1]][[gen]]
      rem1 = table(factor(degree(g),levels = 0 : 50))[-1]
      df1=data.frame(x=as.numeric(names(rev(cumsum(rev(rem1))))), y = rev(cumsum(rev(rem1)))/sum(rem1))
      df1$type = "None removed"
      
      rem2s=replicate(sims,{g2= delete_vertices(g.Ln.f[[fish]][[1]][[gen]],sample(which(MF.aa.f[[fish]][[gen]]$pub >=thr.pub),remove))
      rem2 = table(factor(degree(g2),levels = 0 :50))[-1]})
      df2=data.frame(do.call(rbind,apply(rem2s,2,function(r)data.frame(x=as.numeric(names(rev(cumsum(rev(r))))), y = rev(cumsum(rev(r)))/sum(r)))),curva = rep(1:sims,each = 50))
      df2$type = "Public removed"
      rem3s=replicate(sims,{g3= delete_vertices(g.Ln.f[[fish]][[1]][[gen]],sample(which(MF.aa.f[[fish]][[gen]]$pub <thr.pub),remove))
      rem3 = table(factor(degree(g3),levels = 0 :50))[-1]})
      df3=data.frame(do.call(rbind,apply(rem3s,2,function(r)data.frame(x=as.numeric(names(rev(cumsum(rev(r))))), y = rev(cumsum(rev(r)))/sum(r)))),curva = rep(1:sims,each = 50))
      df3$type = "Non public removed"
      
      deg.maxs[[fish]] = data.frame(deg.max= c(rep(which.min(df1$y[df1$y>0]),40),tapply(df2$y,df2$curva,function(x)which.min(x[x>0])),tapply(df3$y,df3$curva,function(x)which.min(x[x>0]))),cond = factor(rep(c("None removed","Public removed","Non public removed"),each = 40), levels = c("None removed","Public removed","Non public removed")))
      
      p[[fish]]=  ggplot() + labs(x="log(k)", y="log(CDF)") + theme_bw() + 
        geom_line(data = df2, aes(x=log10(x), y=log10(y), group = curva,color= type)) +
        geom_line(data = df3, aes(x=log10(x), y=log10(y), group = curva,color= type)) + 
        geom_line(data = df1, aes(x=log10(x), y=log10(y), group = type, color= type)) + scale_color_manual(values = c("light blue","black","dark red")) +
        ggtitle (label = paste("Removing (2/3 of public clones)",remove,"of", length(V(g)), "nodes ", MF.fish.names[[fish]])) +
        xlab ("log 10 degree") + ylab ("log 10 cumulative degree frequency")
      # dev.off()
    }
  }
  #pdf(file = "Minifish paper/Node removal IgH by fish.pdf",width = 12,height = 8)
  grid.arrange(grobs = p) 
  #dev.off()
  }
  #pdf(file = "Minifish paper/Node removal IgH by fish boxplots.pdf",width = 10,height = 5)
  
  df = data.frame(do.call(rbind,deg.maxs),fish = rep(MF.fish.names,each = 120))
  ggplot (df) + geom_boxplot(aes(y = deg.max,group = c(cond),colour = factor(cond))) + 
    facet_wrap(~fish,nrow = 1) + guides(fill=guide_legend(title="Publicity")) + ylab ("maximum degree") + 
    guides(color=guide_legend(title="Node Removal")) +
    scale_color_manual(values = c("black","dark red","light blue")) +
    theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
  #dev.off()
  
  sapply(1:4,function(i)wilcox.test(data.frame(df%>%filter(cond == "Public removed" & fish == MF.fish.names[[i]]))$deg.max,data.frame(df%>%filter(cond == "Non public removed" & fish == MF.fish.names[[i]]))$deg.max)$p.value)
  
  #write.csv2(df,file = "Minifish paper/Node removal IgH by fish.source.csv")
# Removing all genes for fish 4
  
{ p = list()
  deg.maxs = list()
  for (fish in 4:4){
    for (gen in 1:5){
      thr.pub = 2
      sims = 40
      remove = floor(sum(MF.aa.f[[fish]][[gen]]$pub>=thr.pub)*1)
      #remove = 20
      #remove = 200  
      g = g.Ln.f[[fish]][[1]][[gen]]
      rem1 = table(factor(degree(g),levels = 0 : 50))[-1]
      df1=data.frame(x=as.numeric(names(rev(cumsum(rev(rem1))))), y = rev(cumsum(rev(rem1)))/sum(rem1))
      df1$type = "None removed"
      
      rem2s=replicate(sims,{g2= delete_vertices(g.Ln.f[[fish]][[1]][[gen]],sample(which(MF.aa.f[[fish]][[gen]]$pub >=thr.pub),remove))
      rem2 = table(factor(degree(g2),levels = 0 :50))[-1]})
      df2=data.frame(do.call(rbind,apply(rem2s,2,function(r)data.frame(x=as.numeric(names(rev(cumsum(rev(r))))), y = rev(cumsum(rev(r)))/sum(r)))),curva = rep(1:sims,each = 50))
      df2$type = "Public removed"
      rem3s=replicate(sims,{g3= delete_vertices(g.Ln.f[[fish]][[1]][[gen]],sample(which(MF.aa.f[[fish]][[gen]]$pub <thr.pub),remove))
      rem3 = table(factor(degree(g3),levels = 0 :50))[-1]})
      df3=data.frame(do.call(rbind,apply(rem3s,2,function(r)data.frame(x=as.numeric(names(rev(cumsum(rev(r))))), y = rev(cumsum(rev(r)))/sum(r)))),curva = rep(1:sims,each = 50))
      df3$type = "Non public removed"
      
      deg.maxs[[gen]] = data.frame(deg.max= c(rep(which.min(df1$y[df1$y>0]),40),tapply(df2$y,df2$curva,function(x)which.min(x[x>0])),tapply(df3$y,df3$curva,function(x)which.min(x[x>0]))),cond = factor(rep(c("None removed","Public removed","Non public removed"),each = 40), levels = c("None removed","Public removed","Non public removed")))
      
      p[[gen]]=  ggplot() + labs(x="log(k)", y="log(CDF)") + theme_bw() + 
        geom_line(data = df2, aes(x=log10(x), y=log10(y), group = curva,color= type)) +
        geom_line(data = df3, aes(x=log10(x), y=log10(y), group = curva,color= type)) + 
        geom_line(data = df1, aes(x=log10(x), y=log10(y), group = type, color= type)) + scale_color_manual(values = c("light blue","black","dark red")) +
        ggtitle (label = paste("Removing (all pc)",remove,"of", length(V(g)), "nodes ", MF.gene.names[[gen]])) +
        xlab ("log 10 degree") + ylab ("log 10 cumulative degree frequency")
      # dev.off()
    }
  }
  #pdf(file = "Minifish paper/Node removal by gene.pdf",width = 10,height = 8)
  grid.arrange(grobs = p) 
  #dev.off()
}  
  #pdf(file = "Minifish paper/Node removal by fish boxplots.pdf",width = 12,height = 4)
  
  df = data.frame(do.call(rbind,deg.maxs),gen = factor(rep(MF.gene.names,each = 120),levels = MF.gene.names))
  ggplot (df) + geom_boxplot(aes(y = deg.max,group = c(cond),colour = factor(cond))) + 
    facet_wrap(~gen,nrow = 1) + guides(fill=guide_legend(title="Publicity")) + ylab ("maximum degree") + 
    guides(color=guide_legend(title="Node Removal")) +
    scale_color_manual(values = c("black","dark red","light blue")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
     theme(strip.text = element_text(face = "italic"))
  
  #dev.off()
  #write.csv2(df,file = "Minifish paper/Node removal by fish.source.csv")

# Alpha graph
{  
  lab = NA
  gen=1
  g.a.f5 = g.Ln.f[[4]][[1]][[gen]]
  l.a.f5 = layout.kamada.kawai(g.a.f5)
  js.colors.f5 = gg_color_hue(max(as.numeric(unique(MF.aa.f[[4]][[gen]]$J)))+1)[1+as.numeric(MF.aa.f[[4]][[gen]]$J)]
  js.colors.f5[MF.aa.f[[4]][[gen]]$J ==0] = "white"
  #pdf("Minifish paper//network plots.TCRa.pdf")
  plot(g.a.f5,layout =l.a.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()
}
# Beta graph
{
  fish = 4
  gen = 2
  g.b.f5 = g.Ln.f[[fish]][[1]][[gen]]
  l.b.f5 = layout.kamada.kawai(g.b.f5)
  js.colors.f5 = gg_color_hue(max(as.numeric(unique(MF.aa.f[[fish]][[gen]]$J))))[1+as.numeric(MF.aa.f[[fish]][[gen]]$J)]
  js.colors.f5[MF.aa.f[[fish]][[gen]]$J ==0] = "white"
  #pdf("Minifish paper//network plots.TCRb.pdf")
  plot(g.b.f5,layout =l.b.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()
  
  vs.colors.f5 = gg_color_hue(max(as.numeric(factor(unique(MF.aa.f[[fish]][[gen]]$V)))))[as.numeric(factor(MF.aa.f[[fish]][[gen]]$V))]
  vs.colors.f5[MF.aa.f[[fish]][[gen]]$V ==0] = "white"
  #pdf("Minifish paper//network plots.TCRb.V.pdf")
  plot(g.b.f5,layout =l.b.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()
}
# Gamma graph
  {
  fish = 4
  gen = 3
  g.g.f5 = g.Ln.f[[fish]][[1]][[gen]]
  #to.remove = sample(1:length(V(g.d.f5)),300)
  #g.g.f5 = delete_vertices(g.g.f5,to.remove)
  l.g.f5 = layout.kamada.kawai(g.g.f5)
  js.colors.f5 = gg_color_hue(max(as.numeric(unique(MF.aa.f[[fish]][[gen]]$J))))[as.numeric(MF.aa.f[[fish]][[gen]]$J)]
  js.colors.f5[MF.aa.f[[fish]][[gen]]$J ==0] = "white"
  #pdf("Minifish paper/network plots.TCRg.pdf")
  plot(g.g.f5,layout =l.g.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()
  }
# Delta graph
  {
  fish = 4
  gen = 4
  g.d.f5 = g.Ln.f[[fish]][[1]][[gen]]
  #to.remove = sample(1:length(V(g.d.f5)),300)
  #g.d.f5 = delete_vertices(g.d.f5,to.remove)
  l.d.f5 = layout.kamada.kawai(g.d.f5)
  js.colors.f5 = gg_color_hue(max(as.numeric(unique(MF.aa.f[[fish]][[gen]]$J))))[as.numeric(MF.aa.f[[fish]][[gen]]$J)]
  js.colors.f5[MF.aa.f[[fish]][[gen]]$J ==0] = "white"
  #pdf("Minifish paper//network plots.TCRd.pdf")
  plot(g.d.f5,layout =l.d.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()
  }
# IgH graph
  fish = 4
  gen = 5
  g.h.f5 = g.Ln.f[[fish]][[1]][[gen]]
  l.h.f5 = layout.kamada.kawai(g.h.f5)
  js.colors.f5 = gg_color_hue(max(as.numeric(unique(MF.aa.f[[fish]][[gen]]$J)))+1)[1+as.numeric(MF.aa.f[[fish]][[gen]]$J)]
  js.colors.f5[MF.aa.f[[fish]][[gen]]$J ==0] = "white"
  #pdf("Minifish paper/network plots.IgH.pdf")
  plot(g.h.f5,layout =l.h.f5,vertex.label = lab, vertex.label.color = "black", vertex.label.cex = 1.5, vertex.color = js.colors.f5,vertex.size = MF.aa.f[[4]][[gen]]$pub)
  #dev.off()  
  
# General Statistics table for all Graphs at distance 1  
{
  df.stats = data.frame(gene = factor(rep(MF.gene.names,4),levels = MF.gene.names),fish = rep(MF.fish.names,each=5))
  df.stats$edges = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)length(E(g.Ln.f[[j]][[1]][[i]])))))
  df.stats$nodes = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)length(V(g.Ln.f[[j]][[1]][[i]])))))
  df.stats$mean_degree = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)mean(degree(g.Ln.f[[j]][[1]][[i]])))))
  df.stats$max_degree = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)max(degree(g.Ln.f[[j]][[1]][[i]])))))
  df.stats$max_csize = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)max(clusters(g.Ln.f[[j]][[1]][[i]])$csize))))
  df.stats$mean_csize = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)mean(clusters(g.Ln.f[[j]][[1]][[i]])$csize))))
  
  df.stats$max_diam = as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)diameter(induced_subgraph(g.Ln.f[[j]][[1]][[i]],which(clusters(g.Ln.f[[j]][[1]][[i]])$membership%in%which.max(clusters(g.Ln.f[[j]][[1]][[i]])$csize)))))))
  df.stats$mean_betweenness =as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)mean(MF.aa.f[[j]][[i]]$betweenness))))
  df.stats$mean_authority =as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)mean(MF.aa.f[[j]][[i]]$authority))))
  df.stats$mean_coreness =as.numeric(sapply(1:4,function(j)sapply(1:5,function(i)mean(MF.aa.f[[j]][[i]]$coreness))))
} 

  df.stats
