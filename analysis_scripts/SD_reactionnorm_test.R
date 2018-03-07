setwd("~/raw_count_data/")

library("ggplot2")
library("DESeq2")
library("plyr")
library("dplyr")
library("tidyr")
library("tibble")
library("VennDiagram")
library("seqinr")
library("gridExtra")
library("AnnotationDbi")
library("org.Dm.eg.db")


experimental_design <- read.csv("express_design_matrix.csv",row.names = 1)


readcounts_full <- read.csv("full_index_reads_rounded_annot.csv", row.names = 1)



readcounts_horn <- readcounts_full[, grep(".*hdhorn", colnames(readcounts_full),)]

expdesign_horn <- experimental_design[grep("hdhorn", experimental_design$tissue),]
colnames(readcounts_horn) == row.names(expdesign_horn)


horn.dds <- DESeqDataSetFromMatrix(readcounts_horn, expdesign_horn , 
                                   design = formula(~ sex + size + size:sex))

horn.dds <-DESeq(horn.dds, quiet=TRUE)



##thoracic horn
readcounts_th <- readcounts_full[, grep(".*thxhorn", colnames(readcounts_full),)]

expdesign_th <- experimental_design[grep("thxhorn", experimental_design$tissue),]
colnames(readcounts_th) == row.names(expdesign_th)

th.dds <- DESeqDataSetFromMatrix(readcounts_th, expdesign_th , 
                                 design = formula(~ sex + size + size:sex))

th.dds <-DESeq(th.dds, quiet=TRUE)


##wings
readcounts_win <- readcounts_full[, grep(".*wings", colnames(readcounts_full),)]

expdesign_win <- experimental_design[grep("wings", experimental_design$tissue),]
colnames(readcounts_win) == row.names(expdesign_win)

win.dds <- DESeqDataSetFromMatrix(readcounts_win, expdesign_win , 
                                  design = formula(~ sex + size + size:sex))

win.dds <-DESeq(win.dds, quiet=TRUE)


horn.de<-subset(results(horn.dds, alpha=0.05), padj <0.05)
th.de<-subset(results(th.dds, alpha=0.05), padj <0.05)
win.de<-subset(results(win.dds, alpha=0.05), padj <0.05)


hi.counts <-t(abs(log2((counts(horn.dds[rownames(horn.de), ], 
                               normalized=TRUE, replaced=FALSE)+0.5)))) %>%
  merge(colData(horn.dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(rownames(horn.de))+1):ncol(.)) %>%
  rename(samples=Row.names) %>%
  mutate(samples=as.character(samples))%>%
  mutate(gene=factor(gene))



ti.counts <- t(abs(log2((counts(th.dds[rownames(th.de), ],
                                normalized=TRUE, replaced=FALSE)+0.5)))) %>%
  merge(colData(th.dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(rownames(th.de))+1):ncol(.))%>%
  rename(samples=Row.names) %>%
  mutate(samples=as.character(samples))%>%
  mutate(gene=factor(gene))



wi.counts <- t(abs(log2((counts(win.dds[rownames(win.de), ],
                                normalized=TRUE, replaced=FALSE)+0.5)))) %>%
  merge(colData(win.dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(rownames(win.de))+1):ncol(.))%>%
  rename(samples=Row.names) %>%
  mutate(samples=as.character(samples))%>%
  mutate(gene=factor(gene))




interact.counts <- bind_rows(list(horn = hi.counts, thx = ti.counts, wing=wi.counts), .id = "compare")

summarise.count <-(group_by(interact.counts, gene, sex, size, tissue) %>% 
                     dplyr::summarise(mean=mean(expression)) %>% 
                     tidyr::unite("condition", size , sex, sep=""))

##need to work on this below
ggplot(summarise.count, aes(x=condition, y=mean, color=tissue))+
  geom_point()+geom_line(aes(group=interaction(gene)))+
  #scale_x_discrete(labels=c("Head Horn", "Thorax Horn", "Wings"))+
  labs(x="condition", 
       y="Expression (log normalized counts)")+
  facet_grid(~tissue)




summarise.mean <-interact.counts %>% group_by(gene, sex, size, tissue) %>% 
                    dplyr::summarise(mean=mean(expression)) %>% 
                    unite("condition", size , sex, sep="") %>% 
                    group_by(tissue, condition)%>% 
                    dplyr::summarise(tissuemean=mean(mean)) %>% 
                    mutate(condition = factor(condition, levels=c("lgmale", "smmale", "lgfemale", "smfemale")))

ggplot(summarise.mean, aes(x=condition, y=tissuemean, color=tissue))+
  geom_point()+
  geom_line(aes(group=interaction(tissue)))+
  labs(x="condition", 
       y="Mean Expression (log normalized counts)")





###Standard deviations



all.dds <- DESeqDataSetFromMatrix(readcounts_full, experimental_design , 
                                  design = formula(~ tissue + size + tissue:size))
all.dds$group<- factor(paste0(all.dds$tissue, all.dds$size, all.dds$sex))

design(all.dds) <- ~ group
all.dds <-DESeq(all.dds, quiet=TRUE)

lgvssm.m.hd<-results(all.dds, contrast=c("group", "hdhornlgmale", "hdhornsmmale"), alpha = 0.05)
lgvssm.m.th<-results(all.dds, contrast=c("group", "thxhornlgmale", "thxhornsmmale"), alpha = 0.05)
lgvssm.m.win<-results(all.dds, contrast=c("group", "wingslgmale", "wingssmmale"), alpha = 0.05)
lgvssm.m.gen<-results(all.dds, contrast=c("group", "genitalialgmale", "genitaliasmmale"), alpha = 0.05)


mvf.hd<-results(all.dds, contrast=list(c("grouphdhornlgmale","grouphdhornsmmale"), c("grouphdhornsmfemale","grouphdhornlgfemale")), listValues=c(1/2,-1/2), alpha = 0.05)
mvf.th<-results(all.dds, contrast=list(c("groupthxhornlgmale","groupthxhornsmmale"), c("groupthxhornsmfemale","groupthxhornlgfemale")), listValues=c(1/2,-1/2), alpha = 0.05)
mvf.win<-results(all.dds, contrast=list(c("groupwingslgmale","groupwingssmmale"), c("groupwingssmfemale","groupwingslgfemale")), listValues=c(1/2,-1/2), alpha = 0.05)


size.male.hd.sig <- as.data.frame(subset(lgvssm.m.hd, padj < 0.05))
size.male.th.sig <- as.data.frame(subset(lgvssm.m.th, padj < 0.05))
size.male.win.sig <- as.data.frame(subset(lgvssm.m.win, padj < 0.05))
size.male.gen.sig <- as.data.frame(subset(lgvssm.m.gen, padj < 0.05))


mvf.hd.sig<-as.data.frame(subset(mvf.hd, padj < 0.05))
mvf.th.sig<-as.data.frame(subset(mvf.th, padj < 0.05))
mvf.win.sig<-as.data.frame(subset(mvf.win, padj < 0.05))



de.size<- unique(c(rownames(size.male.hd.sig), rownames(size.male.th.sig), rownames(size.male.win.sig), rownames(size.male.gen.sig)))

de.sex<- unique(c(rownames(mvf.hd.sig), rownames(mvf.th.sig), rownames(mvf.win.sig)))


rld <- rlog(all.dds)
rld.assay <-(assay(rld))


##

de.counts <-t(abs(log2((counts(all.dds[de.size, ], 
                               normalized=TRUE, replaced=FALSE)+0.5)))) %>%
  merge(colData(all.dds), ., by="row.names") %>%
  gather(gene, log2, (ncol(.)-length(de.size)+1):ncol(.)) %>%
  dplyr::rename(samples=Row.names) %>%
  mutate(samples=as.character(samples))%>%
  mutate(gene=factor(gene))



degenes.rld<-t(rld.assay[de.size, ]) %>%
  merge(colData(all.dds), ., by="row.names") %>%
  gather(gene, rld, (ncol(.)-length(de.size)+1):ncol(.)) %>%
  dplyr::rename(samples=Row.names) %>%
  mutate(samples=as.character(samples))%>%
  mutate(gene=factor(gene))


test<- degenes.rld %>% 
  group_by(gene, size, sex) %>% 
  mutate(stdev=sd(rld)) %>% 
  group_by(gene) %>% 
  mutate(mean.sd=mean(stdev)) %>%
  merge(., de.counts)

##this plot didn't take the mean of each log2 for the gene by the two groups
ggplot(test, aes(x=mean.sd, y=log2, shape=size, color=sex))+geom_point(position="jitter")
       

test.2<-degenes.rld %>% 
  group_by(gene, size, sex) %>% 
  mutate(stdev=sd(rld)) %>% 
  group_by(gene) %>% 
  mutate(mean.sd=mean(stdev)) %>% merge(., (de.counts %>% 
                                              group_by(gene, size, sex) %>% 
                                              mutate(mean=mean(log2)) %>% 
                                              group_by(gene) %>% 
                                              mutate(mean.log2=mean(mean))))


#this one takes the mean log2 and mean.sd for each gene across the four "groups" i.e. large males, large females, small males, small females
ggplot(test.2, aes(x=mean.sd, y=mean.log2))+geom_point(position="jitter")

#it looks like the genes with the HIGHEST log2, have the LOWEST SD which I think nicely refutes that comment, unless I misunderstand
 

##finally, testing all of this on the raw counts, instead of combining the RLD object with the raw counts

test.3 <- de.counts %>%
  group_by(gene, size, sex) %>% 
  mutate(mean=mean(log2), sd=sd(log2)) %>% 
  group_by(gene) %>% 
  mutate(mean.log2=mean(mean), mean.sd=mean(sd))

ggplot(test.3, aes(x=mean.sd, y=mean.log2))+geom_point(position="jitter")+geom_smooth(method="lm")
#it looks basically the same so

summary(lm(mean.log2~mean.sd, data=test.3))


#### testing ground


mean.ci<-group_by(interact.counts, gene, sex, size, tissue) %>% 
  dplyr::summarise(mean.ex = mean(expression),
                   sd.ex = sd(expression),
                   n.ex = n()) %>%
  mutate(se.ex = sd.ex / sqrt(n.ex),
         lower.ci = mean.ex - qt(1 - (0.05 / 2), n.ex - 1) * se.ex,
         upper.ci = mean.ex + qt(1 - (0.05 / 2), n.ex - 1) * se.ex)%>%
  unite("condition", size , sex, sep="") %>% 
  group_by(tissue, condition)%>% 
  dplyr::summarise(tissuemean=mean(mean.ex), mean.l.ci=mean(lower.ci), mean.h.ci=mean(upper.ci)) %>% 
  mutate(condition = factor(condition, levels=c("lgmale", "smmale", "lgfemale", "smfemale")))

ggplot(mean.ci, aes(x=condition, y=tissuemean, color=tissue))+
  geom_point()+
  geom_line(aes(group=interaction(tissue)))+
  geom_errorbar(aes(ymin=tissuemean-mean.l.ci, ymax=tissuemean+mean.h.ci), width=.1) +
  geom_point()
labs(x="condition", 
     y="Mean Expression (log normalized counts)")



male.allom<- filter(interact.counts, sex=="male") %>% 
  group_by(individual, tissue) %>%
  dplyr::summarise(mean=mean(expression), mass=log(mean(mass))) 
ggplot(male.allom, aes(x=mass, y=mean, color=tissue))+
  geom_point()+
  geom_smooth()


summary(lm(mean~tissue/mass-1, data=male.allom))
summary(lm(mean~mass*tissue, data=male.allom))
