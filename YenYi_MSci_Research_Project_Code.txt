Project documentation for the large scale analysis of coding and non-coding elements of Mycobacterium tuberculosis, updated as of 11th May 2020.

Step 1: Retrieving Datasets
Datasets are sourced from different databanks such as Array Express, GEOD, and SRA.
For Array Express, simply click on the ENA link at the bottom of the study page. Pressing "Select Columns" option allows the selection of meta data categories, of which we only really need study accession number, run accession number and FASTQ files (FTP). Pressing the "TEXT" option will export the selected sample download links as a .txt file.
For GEOD and SRA, simple extraction is not possible without using in-house tools, which is not easy. However, external software has been created to bypass this problem: https://ewels.github.io/sra-explorer/#. Input the study accession number and the software will automatically query databases and compile a list of samples for the study. Select the desired samples and press "Add to Collection". Once all desired samples are added to collection, press "SAVED DATASETS" (top right corner) and select "Raw FastQ download URLs", which will generate a list of selected sample download links which can be exported as a .txt file.
Move each download link .txt file to the respective study accession directory. To download the files, the command wget is used:
        #!/bin/bash
        wget -i studyaccession.txt
        # This will download the fastq.gz file directly into the directory 
        # Ensure that the .txt file contains the ftp:// prefix
4 datasets are used, including PRJEB65014, PRJNA390669, PRJNA507615 and PRJNA327080. PRJEB19976 was initially used but omitted due to quality issues.    

Step 2: Quality Checks
Quality checks for samples can be done through the software fastQC and multiqc. Multiqc compiles the fastqc results. The modules should be loaded with:
        #!/bin/bash
        # To load fastqc
        module load fastqc

        # To load multiqc with python v3
        module load python/v3

        # To load multiqc with python v2
        module use -s /s/mm/modules
        module load python/v2

The following is a script to run fastQC and multiqc on all files in the directory:
        #!/bin/bash
        # Runs fastqc and multiqc for all .fastq.gz files in a directory
        # Run as:
        # nohup sh fastqc_loop.sh directory_of_samples &

        # Run fastqc for all fastq.gz files in the directory 
        echo "Running fastqc on all files..."
        fastqc *.fastq.gz 

        # Moves output into new folder
        mkdir ./fast_QC_outputs
        mv *fastqc.zip ./fast_QC_outputs
        mv *fastqc.html ./fast_QC_outputs

        # Run multiqc to compile outputs
        cd fast_QC_outputs
        multiqc ./
 
Step 3: Pre-processing
For trimmomatic, the modules should be loaded with:
        #!/bin/bash
        # To load trimmomatic 
        module load trimmomatic

The following is a script to run trimmomatic for paired end samples:
        #!/bin/bash
        # Runs Trimmomatic in PE mode for all sample names given as arguments
        # Run as:
        # nohup sh trimmomatic_PE.sh directory_of_samples list_of_samples

        timestamp=`date "+%Y%m%d-%H%M%S"`
        logfile="run_$timestamp.log"
        exec > $logfile 2>&1  #all output will be logged to logfile

        TRIM_EXEC="/s/software/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar"
        DIR=$1
        shift

        echo "Running Trimmomatic using executable: $TRIM_EXEC"

        for file in `ls $DIR/*_R1_001_trimmed*.fastq.gz` ;
        do
          sample=${file/$DIR\/}
          sample=${sample/_R1_001_trimmed.fastq.gz/}
          echo "Sample= $sample"
          java -jar $TRIM_EXEC PE -threads 12 -phred33 \
               -trimlog "$sample"_trim_report.txt \
               "$DIR$sample"_R1_001_trimmed.fastq.gz "$DIR$sample"_R2_001_trimmed.fastq.gz \
               "$sample"_R1.trimmed_paired.fastq.gz "$sample"_R1.trimmed_unpaired.fastq.gz \
               "$sample"_R2.trimmed_paired.fastq.gz "$sample"_R2.trimmed_unpaired.fastq.gz \
               ILLUMINACLIP:/s/software/trimmomatic/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
               LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

          gzip "$sample"_trim_report.txt
        done
        
Step 4: Mapping of reads to reference genome
This step involves mapping our reads to the reference genome for alignment using the software BWA and samtools. The modules should be loaded with:
        #!/bin/bash

        # To load BWA 
        module load trimmomatic

        # To load samtools
        module load samtools

The reference genome should be indexed with BWA:
        #!/bin/bash
        bwa index ref_genomic.fa

The following is a script to run BWA_mem for paired end samples:
        #!/bin/bash

        # Runs bwa in paired-end  mode 

        # Run as:
        # nohup sh BWA_PE.sh directory_of_fastq_files samples

        timestamp=`date "+%Y%m%d-%H%M%S"`
        logfile="run_$timestamp.log"
        exec > $logfile 2>&1  #all output will be logged to logfile

        dir=$1
        shift

        #set location of executables
        BWA_EXEC="/d/in7/s/bwa-0.7.17/bwa"
        SAMTOOLS_EXEC="/d/in7/s/samtools/samtools-1.3.1/samtools"

        #set parameters
        genomeFile="/d/projects/u/zchayyt/MTB_sample_accession_files/genome_file/ref_genomic.fa
         " #index files should be there too!
        numProc=8

        #extension for fastq files
        suffix1="_R1.trimmed_paired.fastq"
        suffix2="_R2.trimmed_paired.fastq"
        EXT=fastq.gz

        for sample in *.${EXT};
        do
           sample=$(echo $sample | cut -f 1 -d '_')
           echo "Running bwa on sample $sample (paired-end mode)..."

           pairedFile1="$dir$sample$suffix1".gz
           if [ -f $pairedFile1 ]
           then
              gzip -d $pairedFile1
              pairedFile1=$dir$sample$suffix1
           else
              pairedFile1=$dir$sample$suffix1
              if [ ! -f $pairedFile1 ]
              then
                 echo "File not found: $pairedFile1"
                 exit $?
              fi
           fi
           pairedFile2="$dir$sample$suffix2".gz
           if [ -f $pairedFile2 ]
           then
              gzip -d $pairedFile2
              pairedFile2=$dir$sample$suffix2
           else
              pairedFile2=$dir$sample$suffix2
              if [ ! -f $pairedFile2 ]
              then
                 echo "File not found: $pairedFile2"
                 exit $?
              fi
           fi

           tmpSam="$sample"_pe.sam
           tmpBam="$sample"_pe.bam
           finalSortedBam="$sample"_sorted.bam
   
           #align 
           $BWA_EXEC mem -t $numProc $genomeFile $pairedFile1 $pairedFile2 > $tmpSam

           #create bam file
           $SAMTOOLS_EXEC view $tmpSam -Sbo $tmpBam
           $SAMTOOLS_EXEC sort $tmpBam -o $finalSortedBam
           $SAMTOOLS_EXEC index $finalSortedBam

           #cleanup
           /bin/rm $tmpSam $tmpBam
           gzip -9 $pairedFile1 $pairedFile2
        done

Step 5: Mapping output quality check
The following is a script to run mapping output quality checks:
        #!/bin/bash

        timestamp=`date "+%Y%m%d-%H%M%S"`
        logfile="run_$timestamp.log"
        exec > $logfile 2>&1  #all output will be logged to logfile

        dir=$1
        shift

        EXT=bam
        ref_genome="/d/projects/u/zchayyt/MTB_sample_accession_files/genome_file/ref_genomic.bed"
        SUFFIX="_sorted.bam"
        for sample in  *.${EXT};
        do
           sample=$(echo $sample | cut -f 1 -d '_')
           echo "Running mapping quality scripts on sample $sample..."
           echo "sample is $sample"
           quality_check=$dir$sample$SUFFIX
           samtools flagstat $quality_check > "flagstat_$sample.txt"
           echo "Mapping output quality check for $sample done..."
        done

        mkdir flagstat_ouput 
        mv *flagstat* flagstat_output 
        multiqc ./

Step 6: Quantification, normalisation and clustering analysis
The following is an R script used to create a count matrix which undergoes normalisation and hierarchical clustering analysis:
        setwd("~/Desktop/R_Directory")
        #installation of R packages 
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("Rsubread")
        library(Rsubread)
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("DESeq2")
        library(DESeq2)
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("sva")
        library(sva)
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("limma")
        library(limma)
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("WGCNA")
        library(WGCNA)
        install.packages("devtools")
        devtools::install_github("hadley/devtools",force=TRUE)
        library(devtools)
        devtools::install_github("irilenia/baerhunter",force=TRUE)
        library(baerhunter)
        devtools::install_github("zhangyuqing/sva-devel")
        install.packages("ggplot2")
        library(ggplot2)
        install.packages("dendextend")
        library(dendextend)
        install.packages("viridis")
        library(viridis)
        library(colorspace)
        library(RColorBrewer)
        library(Rsamtools)

        #define bam files in the work directory and print file list to confirm samples
        filelist <- list.files(path="./",
                               pattern="*.bam$", 
                               full.names=TRUE, 
                               recursive=FALSE)
        print(filelist)

        #generate count data matrix using featureCounts, for CDS, putative UTRs and putative sRNAs
        countdata_CDS <- featureCounts(files=filelist, 
                                       annot.ext = "./combined_h37rv_Mtb_high_expressed.gff3", 
                                       isGTFAnnotationFile = TRUE, 
                                       chrAliases = "./ChrAlias.txt", 
                                       GTF.attrType = "ID", 
                                       GTF.featureType = "CDS", 
                                       strandSpecific = 2, 
                                       isPairedEnd = TRUE)
        #extracts the counts per feature per sample
        countdata_CDS_final <- countdata_CDS$counts

        countdata_UTR <- featureCounts(files=filelist, 
                                       annot.ext = "./combined_h37rv_Mtb_high_expressed.gff3", 
                                       isGTFAnnotationFile = TRUE, 
                                       chrAliases = "./ChrAlias.txt", 
                                       GTF.attrType = "ID", 
                                       GTF.featureType = "putative_UTR", 
                                       strandSpecific = 2, 
                                       isPairedEnd = TRUE)
        countdata_UTR_final <- countdata_UTR$counts

        countdata_sRNA <- featureCounts(files=filelist, 
                                        annot.ext = "./combined_h37rv_Mtb_high_expressed.gff3", 
                                        isGTFAnnotationFile = TRUE, 
                                        chrAliases = "./ChrAlias.txt", 
                                        GTF.attrType = "ID", 
                                        GTF.featureType = "putative_sRNA", 
                                        strandSpecific = 2, 
                                        isPairedEnd = TRUE)
        countdata_sRNA_final <- countdata_sRNA$counts

        #creates a combined count matrix for all features 
        countdata <- rbind(countdata_CDS_final,countdata_UTR_final,countdata_sRNA_final)
        colnames(countdata) <- gsub("\\.sorted.[sb]am$","",colnames(countdata))
        #inspection of matrix 
        class(countdata)
        head(countdata)
        tail(countdata)

        #defining conditions per sample in a metadata csv file
        control.metadata <- read.csv("control_metadata.csv")
        #if the metadata is somehow already ordered, this line unorders it 
        control.metadata$Condition <- factor(control.metadata$Condition, ordered = FALSE)
        #this step defined control samples as the base reference 
        control.metadata$Condition <- relevel(control.metadata$Condition, "control") 
        control.condition <- factor(control.metadata$Condition)

        #defining coldata for DESeq2
        control.coldata <- data.frame(row.names=colnames(countdata), control.condition)

        #defining dds 
        dds <- DESeqDataSetFromMatrix(countData=countdata, 
                                      colData=control.coldata, 
                                      design=~control.condition)
        #VST transformation; setting vsd parameter of blind = TRUE will result in over-estimating the dispersion and may result in overshriking transformation, resulting in loss of gene signals 
        vsd<-vst(dds,blind=FALSE)

        #boxplot for original data 
        dds.untransformed <- assay(dds)
        colnames(dds.untransformed)<-colnames(dds)
        par(cex.axis=0.5) 
        par(mar=c(4,2,1,1))
        colors = c(rep("#440154FF",3),rep("#31688EFF",22),rep("#FDE725FF",15),rep("#35B779FF",12))
        original.boxplot <- boxplot(dds.untransformed, 
                                    PchCex =0.01,
                                    axes=TRUE,
                                    las=2,
                                    col=colors, 
                                    ylim = c(0,4500),
                                    outline =TRUE,
                                    outcex=0.35)
        legend("topleft", 
               legend=c("E-MTAB-6011", "GEO:GSE67035", "GEO:GSE83814", "GEO:GSE100097"), 
               col = c("#440154FF","#31688EFF","#FDE725FF","#35B779FF"), 
               fill = cols, 
               cex = 4, 
               pt.cex = 1)

        #boxplot for simple logged data
        logged.dds <- log2(countdata + 1)
        colnames(logged.dds)<-colnames(dds)
        par(cex.axis=0.5)
        par(mar=c(4,2,1,1))
        colors = c(rep("#440154FF",3),rep("#31688EFF",22),rep("#FDE725FF",15),rep("#35B779FF",12))
        logged.dds.boxplot <- boxplot(logged.dds,
                                      PchCex =0.01,
                                      axes=TRUE,
                                      las=2,
                                      col=colors,
                                      outcex=0.35)
        legend("topleft", 
               legend=c("PRJEB65014", "PRJNA278760", "PRJNA390669", "PRJNA327080"), 
               col=c("#440154FF","#31688EFF","#FDE725FF","#35B779FF"),
               lty=1, 
               cex=0.75)

        #boxplot for vst data 
        dds.vsd <- assay(vsd)
        colnames(dds.vsd)<-colnames(dds)
        par(cex.axis=0.5)
        par(mar=c(4,2,1,1))
        colors = c(rep("#440154FF",3),rep("#31688EFF",22),rep("#FDE725FF",15),rep("#35B779FF",12))
        vsd.dds.boxplot <- boxplot(dds.vsd,
                                   PchCex =0.01,
                                   axes=TRUE,
                                   las=2,
                                   col=colors,
                                   outcex=0.35)
        legend("topleft", 
               legend=c("PRJEB65014", "PRJNA278760", "PRJNA390669", "PRJNA327080"), 
               col=c("#440154FF","#31688EFF","#FDE725FF","#35B779FF"),
               lty=1, 
               cex=0.75)

        #PCA plotting with DESeq2 in built function 
        PCA.prelim <- plotPCA(vsd,intgroup="control.condition")

        #generating PCA table 
        PCA.data <- data.frame(row.names=colnames(countdata),
                               condition=factor(control.metadata$Condition),
                               dataset=factor(control.metadata$Study))
        PCA.data.plot <- prcomp(t(dds.vsd))
        df.vst <- as.data.frame(PCA.data.plot$x)
        df.vst$group <- PCA.data$condition
        df.vst$dataset <-PCA.data$dataset
        summary(PCA.data.plot)

        #establishing a custom viridis colour palette
        toned_down_pal <- c("#FFBF00","#31688EFF","#35B779FF","#CA0020")

        #PCA plot for VST transformed data 
        custom <- ggplot(df.vst,aes(x=PC1,y=PC2,color=dataset,shape=group)) + scale_shape_manual(values = 0:11) + geom_point(size=3) + xlab("PC1 (42%)") + ylab("PC2 (32%)")
        custom <- custom + scale_color_manual(values = toned_down_pal) + theme_bw()
        custom

        #dendrogram for VST transformed data
        sizeGrWindow(12,9)
        par(cex=0.6)
        par(mar=c(5,6,2,0))
        group <-control.metadata$Study 
        n_group <- length(unique(group)) 
        cols <- toned_down_pal 
        col_group <- cols[group] 
        hc <- hclust(dist(t(dds.vsd)),method="average")
        dend <- as.dendrogram(hc) 
        col_group <- col_group[order.dendrogram(dend)] 
        dend <- dend %>% 
          set("labels_colors", col_group) %>% #change label colors to group
          plot(main = "Dendrogram with samples colored by condition")
        legend("topright", 
               legend = levels(group), 
               fill = cols, 
               cex = 2, 
               pt.cex = 1)

        #batch effect correction using limma; requirement to define batch effect 
        batch.table <- data.frame(study=control.metadata$Study,condition=NA) 
        batch.table$condition <- control.metadata$Condition 
        limma.expr<-removeBatchEffect(x=dds.vsd,
                                      batch=batch.table$study,
                                      batch2=NULL,
                                      covariates=NULL,
                                      design=model.matrix(~ batch.table$condition)) 

        #generating PCA table post-limma
        datExpAdj <- t(limma.expr)
        pca_limma <- prcomp(datExpAdj)
        df.vst.b <- as.data.frame(pca_limma$x)
        df.vst.b$group <- control.metadata$Condition
        df.vst.b$dataset <- control.metadata$Study
        summary(pca_limma) 

        #PCA plot for VST transformed data post-limma; PC1 PC2
        after.limma <- ggplot(df.vst.b,aes(x=PC1,y=PC2,color=dataset, shape=group)) + scale_shape_manual(values = 0:11) + geom_point(size=3) + xlab("PC1 (38%)") + ylab("PC2 (27%)")
        after.limma <- after.limma + scale_color_manual(values = toned_down_pal) + theme_bw()
        after.limma
        #PC2 PC3
        after.limma2 <- ggplot(df.vst.b,aes(x=PC2,y=PC3,color=dataset, shape=group)) + scale_shape_manual(values = 0:11) + geom_point(size=3) + xlab("PC2 (27%)") + ylab("PC3 (22%)")
        after.limma2 <- after.limma2 + scale_color_manual(values = toned_down_pal) + theme_bw()
        after.limma2
        #PC3 PC4
        after.limma3 <- ggplot(df.vst.b,aes(x=PC3,y=PC4,color=dataset, shape=group)) + scale_shape_manual(values = 0:11) + geom_point(size=3) + xlab("PC3 (22%)") + ylab("PC4 (19%)")
        after.limma3 <- after.limma3 + sscale_color_manual(values = toned_down_pal) + theme_bw()
        after.limma3

        #dendrogram for VST transformed data post-limma
        sizeGrWindow(12,9)
        par(cex=0.6)
        par(mar=c(5,6,2,0))
        group <-control.metadata$Study 
        n_group <- length(unique(group)) 
        cols <- toned_down_pal 
        col_group <- cols[group] 
        hc.limma <- hclust(dist(datExpAdj),method="average")
        dend.limma <- as.dendrogram(hc.limma)
        col_group <- col_group[order.dendrogram(dend.limma)] 

        dend.limma <- dend.limma %>% 
          set("labels_colors", col_group) %>% #change label colors to GROUP
          plot(main = "Dendrogram with samples colored by condition")
        legend("center", 
               legend = levels(group), 
               fill = cols, 
               cex = 2,
               pt.cex = 1)

        #WGCNA 
        options(stringsAsFactors = FALSE)
        analysis <- as.data.frame(t(limma.expr))

        #double check data 
        dim(analysis)
        head(analysis)[,1:5]
        #defining trait data 
        traitData <- read.csv("trait_data.csv")
        dim(traitData)
        names(traitData)
        #reformatting data to match
        sample.names <- rownames(analysis)
        traitRows <- match(colnames(dds.vsd),traitData$Sample)
        datTraits <- traitData[traitRows,-1]
        #relabel row names with sample names 
        rownames(datTraits) <- traitData[traitRows,1]
        datTraits[is.na(datTraits)] <- 0 
        collectGarbage()

        #choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        #call the network topology analysis function
        sft = pickSoftThreshold(analysis, 
                                powerVector = powers, 
                                verbose = 5)
        #plotting the output 
        sizeGrWindow(9, 5)
        par(mfrow = c(1,1))
        cex1 = 0.9
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], 
             -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",
             ylab="Scale Free Topology Model Fit,signed R^2",
             type="n",
             main = paste("Scale independence"))
        text(sft$fitIndices[,1], 
             -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,
             cex=cex1,
             col="red")
        #R^2 cut-off line
        abline(h=0.85,col="red")
        #mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], 
             sft$fitIndices[,5],
             xlab="Soft Threshold (power)",
             ylab="Mean Connectivity", 
             type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], 
             sft$fitIndices[,5], 
             labels=powers, 
             cex=cex1,
             col="red")

        #create coexpression network and identify modules
        net = blockwiseModules(analysis, corType = "bicor", networkType = "signed",
                               power = 10, TOMType = "signed", minModuleSize = 20,
                               reassignThreshold = 0, mergeCutHeight = 0.25, deepSplit = 2,
                               numericLabels = TRUE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "testTOM", pamRespectsDendro = FALSE,
                               verbose = 3)
        module.table <- table(net$colors)
        write.csv(table(net$colors), file = "clusteringblocks.csv")
        sizeGrWindow(12, 9)
        #convert labels to colours for plotting
        mergedColors <- labels2colors(net$colors)
        #dendrogram and module colors visual plot 
        plotDendroAndColors(net$dendrograms[[1]], 
                            mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        table(moduleColors)
        write.csv(table(moduleColors), file = "clusteringblocks2.csv")
        MEs = net$MEs
        geneTree = net$dendrograms[[1]]
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "networkConstruction-auto.RData")
        #define numbers of genes and samples
        nGenes = ncol(analysis)
        nSamples = nrow(analysis)
        #recalculate MEs with color labels
        MEs0 = moduleEigengenes(analysis, moduleColors)$eigengenes
        MEs = orderMEs(MEs0)

        #construction of biweight midcorrelation and bonferroni corrected significance of module eigengenes and exprimental conditions
        moduleTraitBicor.data = bicorAndPvalue(MEs, datTraits, maxPOutliers=0.05, robustY = FALSE)
        moduleTraitBicor <- moduleTraitBicor.data$bicor
        moduleTraitBicorPvalue = as.data.frame(moduleTraitBicor.data$p)
        #bonferroni correction applied 
        modNames <- substring(colnames(MEs),3)
        control_bc <- p.adjust(moduleTraitBicorPvalue$Control,method="bonferroni")
        bonferroni.adjusted <- as.data.frame(control_bc,row.names=modNames)
        bonferroni.adjusted$Butyrate <- p.adjust(moduleTraitBicorPvalue$Butyrate,method="bonferroni")
        bonferroni.adjusted$Glucose <- p.adjust(moduleTraitBicorPvalue$Glucose,method="bonferroni")
        bonferroni.adjusted$ButyrateAndGlucose <- p.adjust(moduleTraitBicorPvalue$ButyrateGlucose,method="bonferroni")
        bonferroni.adjusted$HighIron <- p.adjust(moduleTraitBicorPvalue$HighIron,method="bonferroni")
        bonferroni.adjusted$LowIron1Day <- p.adjust(moduleTraitBicorPvalue$LowIron1Day,method="bonferroni")
        bonferroni.adjusted$LowIron1Week <- p.adjust(moduleTraitBicorPvalue$LowIron1Week,method="bonferroni")
        bonferroni.adjusted$TyloxapolAcidic <- p.adjust(moduleTraitBicorPvalue$TyloxapolAcidic,method="bonferroni")
        bonferroni.adjusted$HypoxicLatent <- p.adjust(moduleTraitBicorPvalue$HypoxicLatent,method="bonferroni")
        bonferroni.adjusted$DexStat <- p.adjust(moduleTraitBicorPvalue$DextroseAerobicStationary,method="bonferroni")
        bonferroni.adjusted$CholStat <- p.adjust(moduleTraitBicorPvalue$CholestrolAerobicStationary,method="bonferroni")
        BCmoduleTraitBicorPvalue <- as.matrix(bonferroni.adjusted) 
        #displaying the output in a heatmap 
        textMatrix <- paste(signif(moduleTraitBicor, 2), 
                             "\n(",
                             signif(BCmoduleTraitBicorPvalue, 1), 
                             ")", 
                             sep = "")
        dim(textMatrix) <- dim(moduleTraitBicor)
        sizeGrWindow(10,6)
        par(mar = c(8, 12, 2, 2))
        #heatmap plot 
        labeledHeatmap(Matrix <- moduleTraitBicor,
                       xLabels = colnames(datTraits),
                       yLabels = colnames(MEs),
                       ySymbols = colnames(MEs),
                       colorLabels = TRUE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix2,
                       setStdMargins = FALSE,
                       cex.text = 0.5,
                       cex.lab.x = 0.5,
                       cex.lab.y = 0.5,
                       zlim = c(-1,1),
                       main = paste("Module-trait relationships Bonferroni adjusted p-values"))

        #defining variables containing the experimental condition of datTrait
        hypoxic.latent <- as.data.frame(datTraits$HypoxicLatent)
        names(hypoxic.latent) <- "Hypoxic Latent"
        control <- as.data.frame(datTraits$Control)
        names(control) <- "Control"
        butyrate <- as.data.frame(datTraits$Butyrate)
        names(butyrate) <- "Butyrate"
        glucose <- as.data.frame(datTraits$Glucose)
        names(glucose) <- "Glucose"
        butyrate.glucose <- as.data.frame(datTraits$ButyrateGlucose)
        names(butyrate.glucose) <- "Butyrate and Glucose"
        high.iron <- as.data.frame(datTraits$HighIron)
        names(high.iron) <- "High Iron"
        low.iron.day <- as.data.frame(datTraits$LowIron1Day)
        names(low.iron.day) <- "Low Iron 1 Day"
        low.iron.week <- as.data.frame(datTraits$LowIron1Week)
        names(low.iron.week) <- "Low Iron 1 Week"
        acidic <- as.data.frame(datTraits$TyloxapolAcidic)
        names(acidic) <- "Acidic"
        dex.stat <- as.data.frame(datTraits$DextroseAerobicStationary)
        names(dex.stat) <- "Dextrose Aerobic Stationary"
        chol.stat <- as.data.frame(datTraits$CholestrolAerobicStationary)
        names(chol.stat) <- "Cholestrol Aerobic Stationary"

        #after identifying which modules correlate significantly with interesting experimental conditions, the next step is to examine the genes within these modules and look at the GS and MM; to do so, construction of a data table based on the experimental condition of interest is required

        #generating correlation and significance data table for hypoxia 
        bicor_calc <- as.data.frame(bicorAndPvalue(analysis, MEs))
        geneModuleMembership.b <- bicor_calc[.0:26]
        MMPvalue.b <- bicor_calc[,27:52]
        names(geneModuleMembership.b) = paste("MM", modNames, sep="");
        names(MMPvalue.b) = paste("p.MM", modNames, sep="");
        bicor.gene.trait <-bicorAndPvalue(analysis,hypoxic.latent)
        geneTraitSignificance.b <- as.data.frame(bicor.gene.trait[1])
        GSPvalue.b <- as.data.frame(bicor.gene.trait[2]) 
        names(geneTraitSignificance.b) = paste("GS.", names(hypoxic.latent), sep="");
        names(GSPvalue.b) = paste("p.GS.", names(hypoxic.latent), sep="");

        #using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules, plotting GS against MM 
        module = "brown"
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        sizeGrWindow(7, 7);
        par(mfrow = c(1,1));
        verboseScatterplot(abs(geneModuleMembership.b[moduleGenes, column]),
                           abs(geneTraitSignificance.b[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for hypoxia",
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

        #constructing a data frame holding information for the gene ID, module colours, GS, GS p-value, MM and MM p-values; modules are ordered by significance for weight
        geneInfo0.b = data.frame(gene_ID = names(analysis),
                                 moduleColor = moduleColors,
                                 geneTraitSignificance.b,
                                 GSPvalue.b)
        modOrder.b = order(-abs(bicor(MEs, hypoxic.latent)))
        #adding module membership information in the chosen order
        for (mod in 1:ncol(geneModuleMembership.b)){
          oldNames = names(geneInfo0.b)
          geneInfo0.b = data.frame(geneInfo0.b, geneModuleMembership.b[, modOrder.b[mod]],
                                   MMPvalue.b[, modOrder.b[mod]])
          names(geneInfo0.b) = c(oldNames, paste("MM.", modNames[modOrder.b[mod]], sep=""),
                                 paste("p.MM.", modNames[modOrder.b[mod]], sep=""))
        }
        #order genes in the geneInfo variable by module color, then by geneTraitSignificance
        geneOrder.b= order(geneInfo0.b$moduleColor, -abs(geneInfo0.b$GS.Hypoxic.Latent));
        geneInfo.b = geneInfo0.b[geneOrder.b, ]
        write.csv(geneInfo.b, file = "geneInfo.b.csv")

        #generating correlation and significance data table for control
        bicor_calc <- as.data.frame(bicorAndPvalue(analysis, MEs))
        geneModuleMembership.c <- bicor_calc[.0:26]
        MMPvalue.c <- bicor_calc[,27:52]
        names(geneModuleMembership.c) = paste("MM", modNames, sep="");
        names(MMPvalue.c) = paste("p.MM", modNames, sep="");
        bicor.gene.trait <-bicorAndPvalue(analysis,control)
        geneTraitSignificance.c <- as.data.frame(bicor.gene.trait[1])
        GSPvalue.c <- as.data.frame(bicor.gene.trait[2]) 
        names(geneTraitSignificance.c) = paste("GS.", names(control), sep="");
        names(GSPvalue.c) = paste("p.GS.", names(control), sep="");

        #using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules, plotting GS against MM 
        module = "black"
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        sizeGrWindow(7, 7);
        par(mfrow = c(1,1));
        verboseScatterplot(abs(geneModuleMembership.c[moduleGenes, column]),
                           abs(geneTraitSignificance.c[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for control",
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

        #constructing a data frame holding information for the gene ID, module colours, GS, GS p-value, MM and MM p-values; modules are ordered by significance for hypoxia
        geneInfo0.c = data.frame(gene_ID = names(analysis),
                                 moduleColor = moduleColors,
                                 geneTraitSignificance.c,
                                 GSPvalue.c)
        modOrder.c = order(-abs(bicor(MEs, control)))
        #adding module membership information in the chosen order
        for (mod in 1:ncol(geneModuleMembership.c)){
          oldNames = names(geneInfo0.c)
          geneInfo0.c = data.frame(geneInfo0.c, geneModuleMembership.c[, modOrder.c[mod]],
                                   MMPvalue.c[, modOrder.c[mod]]);
          names(geneInfo0.c) = c(oldNames, paste("MM.", modNames[modOrder.c[mod]], sep=""),
                                 paste("p.MM.", modNames[modOrder.c[mod]], sep=""))
        }
        #order genes in the geneInfo variable by module color, then by geneTraitSignificance
        geneOrder.c= order(geneInfo0.c$moduleColor, -abs(geneInfo0.c$GS.Control));
        geneInfo.c = geneInfo0.c[geneOrder.c, ]
        write.csv(geneInfo.c, file = "geneInfo.c.csv")

        #generating correlation and significance data table for low iron 1 day 
        bicor_calc <- as.data.frame(bicorAndPvalue(analysis, MEs))
        geneModuleMembership.d <- bicor_calc[.0:26]
        MMPvalue.d <- bicor_calc[,27:52]
        names(geneModuleMembership.d) = paste("MM", modNames, sep="");
        names(MMPvalue.d) = paste("p.MM", modNames, sep="");
        bicor.gene.trait <-bicorAndPvalue(analysis,low.iron.day)
        geneTraitSignificance.d <- as.data.frame(bicor.gene.trait[1])
        GSPvalue.d <- as.data.frame(bicor.gene.trait[2]) 
        names(geneTraitSignificance.d) = paste("GS.", names(low.iron.day), sep="");
        names(GSPvalue.d) = paste("p.GS.", names(low.iron.day), sep="");

        #using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules, plotting GS against MM 
        module = "royalblue"
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        sizeGrWindow(7, 7);
        par(mfrow = c(1,1));
        verboseScatterplot(abs(geneModuleMembership.d[moduleGenes, column]),
                           abs(geneTraitSignificance.d[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for low iron 1 day",
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

        #constructing a data frame holding information for the gene ID, module colours, GS, GS p-value, MM and MM p-values; modules are ordered by significance for low iron 1 day
        geneInfo0.d = data.frame(gene_ID = names(analysis),
                                 moduleColor = moduleColors,
                                 geneTraitSignificance.d,
                                 GSPvalue.d)
        modOrder.d = order(-abs(bicor(MEs, low.iron.day)));
        #adding module membership information in the chosen order
        for (mod in 1:ncol(geneModuleMembership.d)){
          oldNames = names(geneInfo0.d)
          geneInfo0.d = data.frame(geneInfo0.d, geneModuleMembership.d[, modOrder.d[mod]],
                                   MMPvalue.d[, modOrder.d[mod]]);
          names(geneInfo0.d) = c(oldNames, paste("MM.", modNames[modOrder.d[mod]], sep=""),
                                 paste("p.MM.", modNames[modOrder.d[mod]], sep=""))
        }
        #order genes in the geneInfo variable by module color, then by geneTraitSignificance
        geneOrder.d= order(geneInfo0.d$moduleColor, -abs(geneInfo0.d$GS.Low.Iron.1.Day));
        geneInfo.d = geneInfo0.d[geneOrder.d, ]
        write.csv(geneInfo.d, file = "geneInfo.d.csv")

        #generating correlation and significance data table for low iron 1 week
        bicor_calc <- as.data.frame(bicorAndPvalue(analysis, MEs))
        geneModuleMembership.e <- bicor_calc[.0:26]
        MMPvalue.e <- bicor_calc[,27:52]
        names(geneModuleMembership.e) = paste("MM", modNames, sep="");
        names(MMPvalue.e) = paste("p.MM", modNames, sep="");
        bicor.gene.trait <-bicorAndPvalue(analysis,low.iron.week)
        geneTraitSignificance.e <- as.data.frame(bicor.gene.trait[1])
        GSPvalue.e <- as.data.frame(bicor.gene.trait[2]) 
        names(geneTraitSignificance.e) = paste("GS.", names(low.iron.week), sep="");
        names(GSPvalue.e) = paste("p.GS.", names(low.iron.week), sep="");

        #using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules, plotting GS against MM 
        module = "salmon"
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        sizeGrWindow(7, 7);
        par(mfrow = c(1,1));
        verboseScatterplot(abs(geneModuleMembership.e[moduleGenes, column]),
                           abs(geneTraitSignificance.e[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for low iron 1 week",
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

        #constructing a data frame holding information for the gene ID, module colours, GS, GS p-value, MM and MM p-values; modules are ordered by significance for low iron 1 week 
        geneInfo0.e = data.frame(gene_ID = names(analysis),
                                 moduleColor = moduleColors,
                                 geneTraitSignificance.e,
                                 GSPvalue.e)
        modOrder.e = order(-abs(bicor(MEs, low.iron.week)));
        #adding module membership information in the chosen order
        for (mod in 1:ncol(geneModuleMembership.e)){
          oldNames = names(geneInfo0.e)
          geneInfo0.e = data.frame(geneInfo0.e, geneModuleMembership.e[, modOrder.e[mod]],
                                   MMPvalue.e[, modOrder.e[mod]]);
          names(geneInfo0.e) = c(oldNames, paste("MM.", modNames[modOrder.e[mod]], sep=""),
                                 paste("p.MM.", modNames[modOrder.e[mod]], sep=""))
        }
        #order genes in the geneInfo variable by module color, then by geneTraitSignificance
        geneOrder.e= order(geneInfo0.e$moduleColor, -abs(geneInfo0.e$GS.Low.Iron.1.Week));
        geneInfo.e = geneInfo0.e[geneOrder.e, ]
        write.csv(geneInfo.e, file = "geneInfo.e.csv")

        ####

