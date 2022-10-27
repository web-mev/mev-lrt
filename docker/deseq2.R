suppressMessages(suppressWarnings(library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)

# Path to the raw/integer counts (tab-delim):
RAW_COUNT_MATRIX<-args[1]

# Path to an annotation/design matrix (tab-delim)
ANNOTATIONS_MATRIX<-args[2]

# The covariate of interest for our LRT. Should be
# one of the column names in the annotation matrix.
ORIG_COVARIATE<-args[3]

OUTPUT_DESEQ_FILE_BASE <- 'deseq2_LRT_results'
OUTPUT_NORMALIZED_COUNTS <- 'deseq2_normalized_counts.tsv'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# Load the annotations:
annotations = read.table(ANNOTATIONS_MATRIX, sep='\t', header=T, row.names=1, stringsAsFactors=T)

# mutate the covariate name so that it matches the column name from the annotation matrix
COVARIATE = make.names(ORIG_COVARIATE)

# check that the desired covariate is actually given in the annotation file. If not, fail out immediately
if (! (COVARIATE %in% colnames(annotations)) ) {
    message(sprintf('The column "%s" was not found in your annotation file. Please check your spelling.', ORIG_COVARIATE))
    quit(status=1)
}

# create a column of the "R-mangled" sample IDs. Since the import of the count matrix will change those names, we need
# to ensure that our annotation file also changes the names.
# Note that we also mangle the rownames. This way when we reorder (and possibly subset)
# the annotation matrix, we maintain annotations as a dataframe. If we do a row-select
# on a single-column dataframe, you end up with a factor, which isn't what we want
# for creationg the DESeq object below.
SAMPLE_ID_COL <- '__sample_name__' 
annotations[SAMPLE_ID_COL] <- make.names(rownames(annotations))
rownames(annotations) <- make.names(rownames(annotations))

# subset to keep only that covariate of interest and rename it for ease:
annotations <- annotations[c(COVARIATE, SAMPLE_ID_COL)]
colnames(annotations) <- c('group', SAMPLE_ID_COL)

# ensure it is a factor (e.g. in case they encode it as integers):
annotations$group <- as.factor(annotations$group)

n_groups <- length(levels(annotations$group))
# check that we actually have more than 1 level for this covariate:
if(n_groups < 2 ) {
    message(sprintf('The covariate "%s" did not have >= 2 levels; we found: %s. To perform a contrast, you need at least two groups to compare.', 
        ORIG_COVARIATE, paste(levels(annotations$group), sep=',')))
    quit(status=1)
}

LEVEL_MAX = 10
if(n_groups > LEVEL_MAX ) {
    message(sprintf('The covariate "%s" had greater than %s levels; the total was %s. 
        Is it possible you are attempting to test the significance 
        of a continuous covariate? To avoid issues with that, we fail 
        jobs that have greater than %s levels.', 
        ORIG_COVARIATE, n_groups, LEVEL_MAX, LEVEL_MAX))
    quit(status=1)
}

# what is the minimum number of samples per level?
NMIN = 2
level_counts <- table(annotations$group) < NMIN
if(any(level_counts)){
    problem_groups <- paste(names(which(level_counts)), collapse=', ')
    message(sprintf('We require %s or more samples for each group. 
        We found the following groups did not meet this requirement: %s', NMIN, problem_groups))
        quit(status=1)
}

# read the raw count matrix, genes as row names. Note that we block R from automatically mangling the names. This way we can preserve
# the original sample names for export:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F, check.names=F)

original_sample_names = colnames(count_data)
modified_sample_names = make.names(original_sample_names)
colname_mapping = data.frame(
    orig_names = original_sample_names,
    row.names=modified_sample_names,
    stringsAsFactors=F)

# now re-assign the column names.
colnames(count_data) <- modified_sample_names

# intersect the names of the samples from the annotation matrix and those from the count matrix:
common_samples = intersect(colnames(count_data), annotations[,SAMPLE_ID_COL])

# if no samples in common, fail out:
if (length(common_samples) == 0){
    message('We could not find any samples in common between your annotation and count matrix. Please check that they are spelled the same.')
    quit(status=1)
}

# Subset the counts and annotations-- they will be in the same order:
annotations <- annotations[common_samples,]
count_data <- count_data[,common_samples]

# run the actual differential expression:
dds <- DESeqDataSetFromMatrix(countData = count_data,
							  colData = annotations,
							  design = ~group)


# wraps the typical DESeq call to catch edge cases where the 
# table of expressions is small and we cannot use the typical
# methods to estimate the mean-dispersion relationship.
run_dge_func <- function(dds){
    tryCatch(
        {
            dds <- DESeq(dds, test='LRT', reduced= ~ 1)
            return (dds)
        },
        error=function(x){
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
            dispersions(dds) <- mcols(dds)$dispGeneEst
            dds <- nbinomLRT(dds)
            return (dds)
        }
    )
}

dds <- run_dge_func(dds)
res <- results(dds, cooksCutoff=F)

# gets text descriptions of each column
col_desc <- mcols(res)$description

# iterate through the column descriptions so we can extract
# which contrast was responsible for the estimate of LFC:
found <- FALSE
i <- 1
while (!found) {
    x = col_desc[i]
    if (startsWith(x, 'log2 fold change')){
        found <- TRUE
        lfc_comparison <- trimws(strsplit(x, ':')[[1]][2]) # looks like "group X vs Y"
        lfc_comparison <- substr(lfc_comparison, 7, nchar(lfc_comparison))
    }
    i <- i+1
}

if (!found){
    message('An unexpected error occurred. An administrator has been notified')
    quit(status=1)
}

column_subset <- c('baseMean', 'log2FoldChange', 'stat', 'pvalue', 'padj')
res <- res[,column_subset]
res <- cbind(Gene=rownames(res), res)

# extract and output the normalized counts:
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)

# map the potentially 'mangled' names back to the original
nc_cols = colnames(nc)
remapped_cols = colname_mapping[nc_cols, 'orig_names']
colnames(nc) = remapped_cols
nc <- cbind(gene=rownames(nc), nc)
fout2 <- paste(working_dir, OUTPUT_NORMALIZED_COUNTS, sep='/')
write.table(nc, fout2, sep='\t', quote=F, row.names=F)

# merge to create a single table, which makes frontend work easier
m <- merge(res, nc, by.x="Gene", by.y=0)
m <- m[order(m$pvalue),]
drops <- c("gene")
m <- m[, !(names(m) %in% drops)]

# change column name for the 'stat' column which will match other dge-type analyses
cols <- colnames(m)
cols[which(names(m) == 'stat')] = 'statistic'
colnames(m) <- cols

output_filename <- paste(OUTPUT_DESEQ_FILE_BASE, COVARIATE, 'tsv', sep='.')
output_filename <- paste(working_dir, output_filename, sep='/')
write.table(m, output_filename, sep='\t', quote=F)

json_str = paste0(
       '{"dge_results":"', output_filename, '",',
       '"lfc_comparison":"', lfc_comparison, '",',
       '"normalized_counts":"', OUTPUT_NORMALIZED_COUNTS, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)

