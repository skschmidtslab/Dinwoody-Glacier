### Analysis of processed amplicon sequence data from Dinwoody Glacier 
  # collected in 2019
setwd("/home/pacifica/R/winds/")
library(tidyverse)
#library(phyloseq)

### Load and format data
  din16s_wtax <- read.table("Din16s_seqtab_wTax_mctoolsr.txt",header = TRUE)
    colnames(din16s_wtax) <- gsub("Undetermined_S0_L001_R1_001.fastq_","",colnames(din16s_wtax))
    tax16s <- data.frame(din16s_wtax$ESV_ID,din16s_wtax$taxonomy)
    rownames(din16s_wtax) <- din16s_wtax$ESV_ID
    analysis16s <- din16s_wtax[,2:31]
    analysis16s <- analysis16s[rowSums(analysis16s) > 0,]
    tax16s_analysis <- tax16s %>% filter(din16s_wtax.ESV_ID %in% rownames(analysis16s))
  din18s_wtax <- read.table("Din18s_seqtab_wTax_mctoolsr.txt",header = TRUE)
    colnames(din18s_wtax) <- gsub("Undetermined_S0_L001_R1_001_","",colnames(din18s_wtax))
    tax18s <- data.frame(din18s_wtax$ESV_ID,din18s_wtax$taxonomy)
    rownames(din18s_wtax) <- din18s_wtax$ESV_ID
    analysis18s <- din18s_wtax[,2:31]
    analysis18s <- analysis18s[rowSums(analysis18s) > 0,]
    tax18s_analysis <- tax18s %>% filter(din18s_wtax.ESV_ID %in% rownames(analysis18s))
  metadata <- read.csv("Metadata_forR.csv") %>% as_tibble()
  
  # Format and filter by taxonomy
  tax16s_wide <- tax16s %>% separate_wider_delim(din16s_wtax.taxonomy, 
                                                 delim = ";", 
                                                 names = c("tax1","tax2","tax3","tax4","tax5","tax6","tax7"))
  tax16s_filt <- tax16s_wide %>% 
    filter(tax1 != "Eukaryota") %>% 
    filter(tax4 != "Chloroplast") %>% 
    filter(tax5 != "Mitochondria") %>% 
    dplyr::select(ESV=tax7,Domain=tax1,Phylum=tax2,Class=tax3,Order=tax4,Family=tax5,Genus=tax6)
  analysis16s <- analysis16s %>% filter(rownames(analysis16s) %in% tax16s_filt$ESV)
  tax16s_filt <- tax16s_filt %>% filter(ESV %in% rownames(analysis16s))
    archaea <- tax16s_filt %>% filter(Domain=="Archaea")
    archaea_count <- analysis16s %>% filter(rownames(analysis16s) %in% archaea$ESV)
    round(100*sum(archaea_count)/sum(analysis16s),2)
  
  tax18s_wide <- tax18s %>% separate_wider_delim(din18s_wtax.taxonomy, 
                                                 delim = ";", 
                                                 names = c("tax1","tax2","tax3","tax4","tax5","tax6","tax7","tax8","tax9"))
  tax18s_filt <- tax18s_wide %>% 
    filter(tax1 != "Bacteria") %>% 
    filter(tax7 != "Homo_sapiens_(human)") %>% 
    dplyr::select(ESV=din18s_wtax.ESV_ID,Domain=tax1,Phylum=tax2,Class=tax3,Order=tax4,Family=tax5,Genus=tax6)
  analysis18s <- analysis18s %>% filter(rownames(analysis18s) %in% tax18s_filt$ESV)
  tax18s_filt <- tax18s_filt %>% filter(ESV %in% rownames(analysis18s))
  

### Compare alpha diversity and plot  
  # Using Breakaway (what's the input? see vignette at https://adw96.github.io/breakaway/articles/breakaway.html)
library(breakaway)
    # test and plot for 16S (bacteria)
    richness_din16s <- breakaway(analysis16s)
      summary(richness_din16s) %>% as_tibble
    combined_richness_16s <- metadata %>% dplyr::select(sample_names=sample,elevation) %>% 
      left_join(summary(richness_din16s))
    bt_elev_16s <- betta(formula = estimate ~ elevation,
                     ses = error, data = combined_richness_16s)
    bt_elev_16s$table
    richnessplot_16s <- ggplot(combined_richness_16s,aes(x=elevation,y=estimate)) + 
                        geom_boxplot(outlier.color = 'white') +
                        geom_point(position=position_dodge(w=.8)) +
                        labs(x="",y="Estimated ESVs",
                             title="a Bacteria") +
                        theme(
                          plot.title = element_text(face="bold",size=12) 
                          ,plot.background = element_blank()
                          ,panel.background = element_blank()
                          ,panel.grid.major = element_blank()
                          ,panel.grid.minor = element_blank()
                          ,panel.border = element_blank()
                          ,axis.line = element_line(colour="black", size=1)
                          ,axis.text.x = element_text(color="black", size=12)
                          ,axis.text.y = element_text(color="black", size=12)
                          ,axis.title.x = element_text(colour="black", size=12)
                          ,axis.title.y = element_text(color="black", size=12)
                          ,legend.position="none"
                        )
    # test and plot for 18S (eukaryotes)
    richness_din18s <- breakaway(analysis18s)
    summary(richness_din18s) %>% as_tibble
    combined_richness_18s <- metadata %>% dplyr::select(sample_names=sample,elevation) %>% 
      left_join(summary(richness_din18s))
    bt_elev_18s <- betta(formula = estimate ~ elevation,
                         ses = error, data = combined_richness_18s)
    bt_elev_18s$table
    richnessplot_18s <- ggplot(combined_richness_18s,aes(x=elevation,y=estimate)) + 
      geom_boxplot(outlier.color = 'white') +
      geom_point(position=position_dodge(w=.8)) +
      labs(x="",y="Estimated ESVs",
           title="b Eukaryotes") +
      theme(
        plot.title = element_text(face="bold",size=12) 
        ,plot.background = element_blank()
        ,panel.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(colour="black", size=1)
        ,axis.text.x = element_text(color="black", size=12)
        ,axis.text.y = element_text(color="black", size=12)
        ,axis.title.x = element_text(colour="black", size=12)
        ,axis.title.y = element_text(color="black", size=12)
        ,legend.position="none"
      )
    # combine and export figure
    library(gridExtra)
    grid.arrange(richnessplot_16s,richnessplot_18s, nrow = 1)
    

### Ordinate and plot community composition
  # Alison used Aitchison ILR and PCA in vegan
    # make less sparse
    # claculate ILR
    # ordinate
    # plot
    # test for significance 
    library(vegan)
    library(zCompositions)
    library(compositions)
    
    # Set and apply our cut-off for mean min reads for keeping taxa:
    field.bac.count <- 20
    field.euk.count <- 20
    analysis16s.20 <- data.frame(analysis16s[which(apply(analysis16s, 1, function(x){mean(x)}) > field.bac.count),], 
                                 check.names=F) #148 ESVs
    analysis18s.20 <- data.frame(analysis18s[which(apply(analysis18s, 1, function(x){mean(x)}) > field.euk.count),], 
                                 check.names=F) #93 ESVs
    analysis16s.czm <- cmultRepl(t(analysis16s.20),  label=0, method="CZM")
    analysis18s.czm <- cmultRepl(t(analysis18s.20),  label=0, method="CZM")
    # REMEMBER: IF samples are columns, margin = 2, if samples are ROWS, margin = 1
    analysis16s.ilr <- ilr(analysis16s.czm)
    analysis18s.ilr <- ilr(analysis18s.czm)
    
    # Make PCA plots
      #Run PCA
      analysis16s.pca <- prcomp(analysis16s.ilr,scale = TRUE)
      analysis18s.pca <- prcomp(analysis18s.ilr,scale = TRUE)
      #Scree plot
      plot(analysis16s.pca,type="l",main="16S (bacteria)")
      plot(analysis18s.pca,type="l",main="18S (euks), count=")
      # Make a readable biplot without the loadings:
        # Sum the total variance
        analysis16s.mvar <- sum(analysis16s.pca$sdev^2)
        analysis18s.mvar <- sum(analysis18s.pca$sdev^2)
        # Calculate the PC1 and PC2 variance
        analysis16s.pc1 <- paste("PC1: ", round(100*sum(analysis16s.pca$sdev[1]^2)/analysis16s.mvar, 1),"%")
        analysis16s.pc2 <- paste("PC2: ", round(100*sum(analysis16s.pca$sdev[2]^2)/analysis16s.mvar, 1),"%")
        analysis18s.pc1 <- paste("PC1: ", round(100*sum(analysis18s.pca$sdev[1]^2)/analysis18s.mvar, 1),"%")
        analysis18s.pc2 <- paste("PC2: ", round(100*sum(analysis18s.pca$sdev[2]^2)/analysis18s.mvar, 1),"%")

      ## 16S plot  
      # Make the blank plot
        plot(analysis16s.pca$x[,2] ~ analysis16s.pca$x[,1], 
             xlab=analysis16s.pc1, 
             ylab=analysis16s.pc2,
             cex=0.1, col="white"
        )
        title("a  Bacteria", adj = 0, line = 0.2)
    # Add the points & groupings
        lowerglacier <- which(metadata$elevation == "Lower glacier")
        upperglacier <- which(metadata$elevation == "Upper glacier")
        Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
          hpts <- chull(x = xcoord, y = ycoord)
          hpts <- c(hpts, hpts[1])
          lines(xcoord[hpts], ycoord[hpts], col = lcolor)
        }
        points(analysis16s.pca$x[lowerglacier,2] ~ analysis16s.pca$x[lowerglacier,1], pch=19,
               cex=1, col="orange")
        Plot_ConvexHull(xcoord = analysis16s.pca$x[lowerglacier,1], ycoord = analysis16s.pca$x[lowerglacier,2], 
                        lcolor = "orange")
        points(analysis16s.pca$x[upperglacier,2] ~ analysis16s.pca$x[upperglacier,1], pch=19,
               cex=1, col="red")
        Plot_ConvexHull(xcoord = analysis16s.pca$x[upperglacier,1], ycoord = analysis16s.pca$x[upperglacier,2], 
                        lcolor = "red")
        
        pca_fp_16s <- tibble(x=analysis18s.pca$x[,2],y=analysis18s.pca$x[,1],Location=metadata$elevation,Sample=metadata$sample)
        pca_16s <- ggplot(pca_fp_16s,aes(x=x,y=y,color=Location)) + 
          geom_point(aes(color=Location)) +
          labs(x="",y="",
               title="a Bacteria") +
          theme(
            plot.title = element_text(face="bold",size=12) 
            ,plot.background = element_blank()
            ,panel.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,axis.line = element_line(colour="black", size=1)
            ,axis.text.x = element_text(color="black", size=12)
            ,axis.text.y = element_text(color="black", size=12)
            ,axis.title.x = element_text(colour="black", size=12)
            ,axis.title.y = element_text(color="black", size=12)
            ,legend.position="none"
          )
        
        
    # 18s plot
        # Make the blank plot
        plot(analysis18s.pca$x[,2] ~ analysis18s.pca$x[,1], 
             xlab=analysis18s.pc1, 
             ylab=analysis18s.pc2,
             cex=0.1, col="white"
        )
        title("b  Eukaryotes", adj = 0, line = 0.2)
        # Add the points & groupings
        points(analysis18s.pca$x[lowerglacier,2] ~ analysis18s.pca$x[lowerglacier,1], pch=19,
               cex=1, col="orange")
        Plot_ConvexHull(xcoord = analysis18s.pca$x[lowerglacier,1], ycoord = analysis18s.pca$x[lowerglacier,2], 
                        lcolor = "orange")
        points(analysis18s.pca$x[upperglacier,2] ~ analysis18s.pca$x[upperglacier,1], pch=19,
               cex=1, col="red")
        Plot_ConvexHull(xcoord = analysis18s.pca$x[upperglacier,1], ycoord = analysis18s.pca$x[upperglacier,2], 
                        lcolor = "red")


        # Then apply a statistical test
        ai.16s <- vegdist(analysis16s.ilr, method = "euclidean",binary = FALSE, diag = FALSE, upper = TRUE)
        adonis2(ai.16s~metadata$elevation) #p=.001 
        bd.16s <- betadisper(ai.16s,metadata$elevation, type = "centroid")
        anova(bd.16s) # p = 0.005, showing difs in dispersion
        
        ai.18s <- vegdist(analysis18s.ilr, method = "euclidean",binary = FALSE, diag = FALSE, upper = TRUE)
        adonis2(ai.18s~metadata$elevation) #p=.001 
        bd.18s <- betadisper(ai.18s,metadata$elevation, type = "centroid")
        anova(bd.18s) # p = 0.77, similar dispersion
        
        
### Plot rank abundance

        # Subset different sampling groups
        lower.bacs <- analysis16s[,1:9]
        upper.bacs <- analysis16s[,10:30]
        lower.euks <- analysis18s[,1:9]
        upper.euks <- analysis18s[,10:30]
        
        # Select top ten taxa by percent abundance for each sampling group & format scientific names
        lower.bacs.per <- 100*sweep(lower.bacs,MARGIN=2,FUN="/",STATS=colSums(lower.bacs))
        lower.bacs.per <- lower.bacs.per[order(rowSums(lower.bacs.per),decreasing = TRUE),]
        lower.bacs.per.10 <- lower.bacs.per[1:10,]
        lower.bacs.per.10$ESV <- rownames(lower.bacs.per.10)
        lower.bacs.per.10.tax <- lower.bacs.per.10 %>% 
          left_join(tax16s_filt) %>% 
          mutate(Genus=gsub("_CCAP_1459-11B","",Genus)) %>% 
          mutate(Genus=gsub("_PCC-7430","",Genus)) %>% 
          mutate(Name = ifelse(Genus=="NA",Family,Genus)) %>%
          mutate(Taxon=paste0(Name," sp. (",ESV,")"))

        upper.bacs.per <- 100*sweep(upper.bacs,MARGIN=2,FUN="/",STATS=colSums(upper.bacs))
        upper.bacs.per <- upper.bacs.per[order(rowSums(upper.bacs.per),decreasing = TRUE),]
        upper.bacs.per.10 <- upper.bacs.per[1:10,]
        upper.bacs.per.10$ESV <- rownames(upper.bacs.per.10)
        upper.bacs.per.10.tax <- upper.bacs.per.10 %>% 
          left_join(tax16s_filt) %>% 
          mutate(Genus=gsub("_CCAP_1459-11B","",Genus)) %>% 
          mutate(Genus=gsub("_PCC-7430","",Genus)) %>% 
          mutate(Name = ifelse(Genus=="NA",Family,Genus)) %>%
          mutate(Name = ifelse(Name=="env.OPS_17",Order,Name)) %>% 
          mutate(Taxon=paste0(Name," sp. (",ESV,")"))
        
        lower.euks.per <- 100*sweep(lower.euks,MARGIN=2,FUN="/",STATS=colSums(lower.euks))
        lower.euks.per <- lower.euks.per[order(rowSums(lower.euks.per),decreasing = TRUE),]
        lower.euks.per.10 <- lower.euks.per[1:10,]
        lower.euks.per.10$ESV <- rownames(lower.euks.per.10)
        lower.euks.per.10.tax <- lower.euks.per.10 %>% 
          left_join(tax18s_filt) %>% 
          mutate(Name = ifelse(Genus=="NA" | Genus=="uncultured" | Genus=="uncultured_eukaryote",Family,Genus)) %>%
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Order,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Class,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Phylum,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote","Unknown eukaryote",Name)) %>% 
          mutate(Taxon=paste0(Name," sp. (",ESV,")"))
        
        upper.euks.per <- 100*sweep(upper.euks,MARGIN=2,FUN="/",STATS=colSums(upper.euks))
        upper.euks.per <- upper.euks.per[order(rowSums(upper.euks.per),decreasing = TRUE),]
        upper.euks.per.10 <- upper.euks.per[1:10,]
        upper.euks.per.10$ESV <- rownames(upper.euks.per.10)
        upper.euks.per.10.tax <- upper.euks.per.10 %>% 
          left_join(tax18s_filt) %>% 
          mutate(Name = ifelse(Genus=="NA" | Genus=="uncultured" | Genus=="uncultured_eukaryote",Family,Genus)) %>%
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Order,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Class,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote",Phylum,Name)) %>% 
          mutate(Name = ifelse(Name=="NA" | Name=="uncultured" | Name=="uncultured_eukaryote","Unknown eukaryote",Name)) %>% 
          mutate(Taxon=paste0(Name," sp. (",ESV,")"))
        
        # Format for plotting and make and combine graphs
        library(reshape2)
        lower.bacs.fp <- melt(lower.bacs.per.10.tax)   
          lower.bacs.fp$Taxon <- reorder(lower.bacs.fp$Taxon, -lower.bacs.fp$value)
        upper.bacs.fp <- melt(upper.bacs.per.10.tax)        
          upper.bacs.fp$Taxon <- reorder(upper.bacs.fp$Taxon, -upper.bacs.fp$value)
        lower.euks.fp <- melt(lower.euks.per.10.tax)        
          lower.euks.fp$Taxon <- reorder(lower.euks.fp$Taxon, -lower.euks.fp$value)
        upper.euks.fp <- melt(upper.euks.per.10.tax)    
          upper.euks.fp$Taxon <- reorder(upper.euks.fp$Taxon, -upper.euks.fp$value)
        
        lower.bacs.plot <- ggplot(lower.bacs.fp, aes(x = Taxon, y = value))  +
          scale_y_continuous(expand = c(0, 0)) +
          geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
          geom_point(aes(color=Phylum)) +
         # scale_fill_manual(c("orange","forestgreen","orange")) + 
          scale_color_manual(values=c("blue","forestgreen","purple")) + 
          #scale_color_gradient(low = "blue", high = "red") +
          labs(x="",y="Percent Abundance",
               title="a   Bacteria on lower glacier") +
          theme(
            plot.title = element_text(color="black",size=12,face="bold") 
            ,plot.background = element_blank()
            ,panel.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,axis.line = element_line(colour="black", size=1)
            ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
            ,axis.text.y = element_text(color="black", size=12)
            ,axis.title.x = element_text(colour="black", size=12)
            ,axis.title.y = element_text(color="black", size=12)
            ,legend.position = "none"
          )        
        lower.bacs.plot
        
        upper.bacs.plot <- ggplot(upper.bacs.fp, aes(x = Taxon, y = value))  +
          #scale_y_continuous(expand = c(0, 0),limits = c(0,21)) +
          scale_color_manual(values=c("blue","forestgreen","purple")) + 
          geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
          geom_point(aes(color=Phylum)) +
          # scale_color_discrete(c("red","forestgreen","red")) +
          labs(x="",y="Percent Abundance",
               title="b   Bacteria on upper glacier") +
          theme(
            plot.title = element_text(color="black",size=12,face="bold") 
            ,plot.background = element_blank()
            ,panel.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,axis.line = element_line(colour="black", size=1)
            ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
            ,axis.text.y = element_text(color="black", size=12)
            ,axis.title.x = element_text(colour="black", size=12)
            ,axis.title.y = element_text(color="black", size=12)
            ,legend.position = "none"
          )        
        upper.bacs.plot
        
        lower.euks.plot <- ggplot(lower.euks.fp, aes(x = Taxon, y = value))  +
          scale_y_continuous(expand = c(0, 0)) +
          geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
          geom_point(aes(color=Phylum)) +
          scale_color_manual(values=c("springgreen4","gray60","maroon4","navy")) + 
          #scale_color_gradient(low = "blue", high = "red") +
          labs(x="",y="Percent Abundance",
               title="c   Eukaryotes on lower glacier") +
          theme(
            plot.title = element_text(color="black",size=12,face="bold") 
            ,plot.background = element_blank()
            ,panel.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,axis.line = element_line(colour="black", size=1)
            ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
            ,axis.text.y = element_text(color="black", size=12)
            ,axis.title.x = element_text(colour="black", size=12)
            ,axis.title.y = element_text(color="black", size=12)
            ,legend.position = "none"
          )        
        lower.euks.plot
        
        upper.euks.plot <- ggplot(upper.euks.fp, aes(x = Taxon, y = value))  +
          scale_y_continuous(expand = c(0, 0)) +
          geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
          geom_point(aes(color=Phylum)) +
          scale_color_manual(values=c("springgreen4","gray60","maroon4","navy")) + 
          labs(x="",y="Percent Abundance",
               title="d   Eukaryotes on upper glacier") +
          theme(
            plot.title = element_text(color="black",size=12,face = "bold") 
            ,plot.background = element_blank()
            ,panel.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,axis.line = element_line(colour="black", size=1)
            ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
            ,axis.text.y = element_text(color="black", size=12)
            ,axis.title.x = element_text(colour="black", size=12)
            ,axis.title.y = element_text(color="black", size=12)
            ,legend.position = "none"
          )        
        upper.euks.plot
        
        ### STILL NEED TO REORDER TAXA IN DECREASING AVERAGE ABUNDANCE!!!!
        grid.arrange(lower.bacs.plot,
                     upper.bacs.plot,
                     lower.euks.plot,
                     upper.euks.plot, 
                     nrow = 2)
        
        
        