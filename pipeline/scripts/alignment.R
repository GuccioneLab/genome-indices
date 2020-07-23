#################################################################
#################################################################
############### Genome Alignment - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. Aggregate counts
#######################################################
#######################################################

#############################################
########## 1. Function to aggregate counts
#############################################
# Aggregates counts from STAR, kallisto or salmon.

aggregate_counts <- function(infiles, method, summarize_to_gene=FALSE, filter_protein_coding=FALSE, name_func=function(x) { x_split <- strsplit(x, '/')[[1]]; x_name <- x_split[length(x_split)-1]; return(x_name); }, star_col=2) {

    # Add names
    names(infiles) <- sapply(infiles, name_func)

    # STAR
    if (method == 'star') {
         
        # Read and concatenate
        concatenated_dataframe <- do.call('rbind', sapply(names(infiles), function(x) {df <- as.data.frame(data.table::fread(infiles[x])); df <- df[grepl('ENSG*', df$V1),]; colnames(df)[1] <- 'ensembl_gene_id'; colnames(df)[star_col] <- 'counts'; df$sample <- x; df}, simplify=FALSE))
        
        # Cast
        result <- reshape2::dcast(ensembl_gene_id ~ sample, data = concatenated_dataframe, value.var='counts', fun.aggregate=sum)

        # Fix rownames
        rownames(result) <- result$ensembl_gene_id
        result$ensembl_gene_id <- NULL
               
    # Kallisto and Salmon
    } else if (method %in% c('kallisto', 'salmon')) {

        # Load
        require(tximeta)
      
        # Create coldata
        coldata <- data.frame(row.names = sapply(infiles, name_func), files = infiles, names = infiles)

        # Get SummarizedExperiment
        result <- tximeta(coldata, type=method)

        # Summarize to gene
        if (summarize_to_gene) {
            result <- summarizeToGene(result)
        }

        # Protein coding
        if (filter_protein_coding) {
            result <- result[rowRanges(result)$gene_biotype == 'protein_coding',]
        }

    }
    
    # Return
    return(result)
        
}

#############################################
########## 2. Tximeta
#############################################
# Aggregates counts from STAR, kallisto or salmon.

import_counts <- function(sample_dataframe, type='salmon', summarize_to_gene=TRUE, filter_protein_coding=TRUE) {

    # Load
    require(tximeta)

    # Get SummarizedExperiment
    se <- tximeta(sample_dataframe, type=type)

    # # Protein coding
    # if (filter_protein_coding) {
    #     se <- se[rowRanges(se)$tx_biotype == 'protein_coding',]
    # }

    # Summarize to gene
    if (summarize_to_gene) {
        se <- summarizeToGene(se)
    }

    # Protein coding
    if (filter_protein_coding) {
        se <- se[rowRanges(se)$gene_biotype == 'protein_coding',]
    }

    # Return
    return(se)
    
}

#############################################
########## 3. Convert to gene symbol
#############################################
# Aggregates and converts counts to gene symbols

convert_gene_symbols <- function(expression_dataframe, genome, merge_col = 'Gene stable ID', protein_coding = FALSE) {

    # Load
    suppressPackageStartupMessages(require(dplyr))

    # Get tx2g
    biomart_dataframe <- as.data.frame(data.table::fread(Sys.glob(glue::glue('/sc/hydra/projects/GuccioneLab/genome-indices/{genome}/ensembl/*-biomart.txt'))))

    # Filter columns
    gene_info_dataframe <- unique(biomart_dataframe[,c(merge_col, 'Gene name', 'Gene type')])

    # Filter PCG
    if (protein_coding) {
        gene_info_dataframe <- gene_info_dataframe[gene_info_dataframe[,'Gene type'] == 'protein_coding',]
    }

    # Merge
    merged_dataframe <- merge(gene_info_dataframe, expression_dataframe, by.x=merge_col, by.y='row.names')

    # Fix columns
    merged_dataframe[,c(merge_col, 'Gene type')] <- NULL
    colnames(merged_dataframe)[1] <- 'gene_symbol'

    # Aggregate
    aggregated_dataframe <- merged_dataframe %>% group_by(gene_symbol) %>% summarize_all(sum) %>% as.data.frame()

    # Fix rownames
    rownames(aggregated_dataframe) <- aggregated_dataframe$gene_symbol
    aggregated_dataframe$gene_symbol <- NULL

    # Return
    return(aggregated_dataframe)
    
}


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################