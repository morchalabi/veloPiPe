library(parallel)
library(Matrix)

# STEP 1: Setting up options ####

cpu_ = 50      # number of CPU cores

# STEP 2: Loading exon and gene tables ####

load(file = '../data/genes_exons_tables.RData')

# STEP 3: Processing batch BAM files

btchs_ = list.files(path = '../data', pattern = 'batch_*', full.names = T)      # reading in batch folders
for(btch_ in btchs_)
{
  message('<< Processing bam files in ',btch_, ' >>')
  
  bamfiles_ = list.files(path = btch_, pattern = '*.bam\\b')                                          # listing bam and bai files
  bamfiles_ = bamfiles_[!grepl(x = bamfiles_, pattern = '*.bai\\b')]                                  # filtering out bai files
  cells_ = sapply(X = strsplit(bamfiles_, split = '[.]'), FUN = function(elm_) return(elm_[[1]]))     # extracting cell names
  bamfiles_ = paste0(btch_,'/',bamfiles_)
  
  rslt_btch = mclapply(X = bamfiles_, mc.cores = cpu_, genes_tbl, exons_tbl,  FUN =
                                                                              function(bamfile_, genes_tbl, exons_tbl)
                                                                              {
                                                                                # step 3.1: Generating alignment table ####
                                                                                
                                                                                tmp_ = paste0(bamfile_,'.txt')                                                                                            # temp out file
                                                                                system2(command = 'samtools', args = c('view', '-b', bamfile_, '|', 'bamToBed', '-i', 'stdin','>', tmp_))
                                                                                aligns_ = read.table(file = tmp_, sep = '\t', header = F, col.names = c('seq','start','end','name','score','strand'))     # alignments are already position sorted
                                                                                system2(command = 'rm', args = tmp_)
                                                                                
                                                                                # step 3.2: Finding spliced and unspliced alignments ####
                                                                                
                                                                                spliced_vec = unspliced_vec = ambiguous_vec = integer(length = nrow(genes_tbl))             # vector containing spliced counts
                                                                                names(spliced_vec) = names(unspliced_vec) = names(ambiguous_vec) =  rownames(genes_tbl)     # vector containing unspliced counts
                                                                                
                                                                                for(seq_ in unique(aligns_$seq))                                    # for each sequence
                                                                                {
                                                                                  genes_tbl_sub = genes_tbl[genes_tbl$seq %in% seq_,]     # table of genes at chr seq_
                                                                                  exons_tbl_sub = exons_tbl[exons_tbl$seq %in% seq_,]     # table of exons at chr seq_
                                                                                  aligns_sub    = aligns_[    aligns_$seq %in% seq_,]     # table of alignments at chr seq_
                                                                                  
                                                                                  for(rw_ in 1:nrow(aligns_sub))                          # for each alignment
                                                                                  {
                                                                                    align_ = aligns_sub[rw_,]                                                                                       # locus of current alignment
                                                                                    align_str_exons = exons_tbl_sub[ exons_tbl_sub$start <= align_$start & align_$start <= exons_tbl_sub$end, ]     # finding exons crossing start of alignment
                                                                                    align_end_exons = exons_tbl_sub[ exons_tbl_sub$start <= align_$end   & align_$end   <= exons_tbl_sub$end, ]     # finding exons crossing end of alignment
                                                                                    genes_ = union(align_str_exons$gene, align_end_exons$gene)                                                      # genes whose exons are crossed by either end of current alignment
                                                                                    
                                                                                    # << Case 1 >>: if neither end of current alignment crosses any exon (alignment is either in intronic or intergenic region)
                                                                                    
                                                                                    if(nrow(align_str_exons) == 0 & nrow(align_end_exons) == 0)
                                                                                    {
                                                                                      gene_ = rownames(genes_tbl_sub[ genes_tbl_sub$start <= align_$start & align_$end <= genes_tbl_sub$end, ])     # gene that contains current alignment
                                                                                      if(length(gene_) == 1){ ambiguous_vec[gene_] = ambiguous_vec[gene_] + 1 }                                     # alignment is in intronic region of a gene
                                                                                      next()
                                                                                    }
                                                                                    
                                                                                    # << Case 2 >>: if both ends of current alignment cross any exon
                                                                                    
                                                                                    if(nrow(align_str_exons) != 0 & nrow(align_end_exons) != 0)
                                                                                    {
                                                                                      cmn_exons = intersect(x = align_str_exons$exon, y = align_end_exons$exon)     # common exons crossed by both ends of current alignment
                                                                                      if(length(cmn_exons) != 0)                                                    # there is at least one exon crossed by both ends of current alignment (alignment is exclusively assigned to such exons)
                                                                                      {
                                                                                        unsp_genes = unique(align_str_exons[align_str_exons$exon %in% cmn_exons,"gene"])      # only genes of common exons
                                                                                        unspliced_vec[unsp_genes] = unspliced_vec[unsp_genes] + 1
                                                                                        
                                                                                      }else                                                                         # crossed exons are apart; they are from either one gene or different genes
                                                                                      {
                                                                                        cmn_genes = intersect(x = align_str_exons$gene, y = align_end_exons$gene)
                                                                                        if(length(cmn_genes) != 0)
                                                                                        {
                                                                                          spliced_vec[cmn_genes] = spliced_vec[cmn_genes] + 1
                                                                                          genes_ = genes_[!genes_ %in% cmn_genes]
                                                                                        }
                                                                                        if(length(genes_) != 0) { ambiguous_vec[genes_] = ambiguous_vec[genes_] + 1 }
                                                                                      }
                                                                                      next()
                                                                                    }
                                                                                    
                                                                                    # << Case 3 >>: only one end of current alignment crosses exon(s)
                                                                                    
                                                                                    ambiguous_vec[genes_] = ambiguous_vec[genes_] + 1
                                                                                  }
                                                                                }
                                                                                
                                                                                system2(command = 'echo', args = c('done for ', bamfile_))
                                                                                return( cbind(unspliced = unspliced_vec, spliced = spliced_vec, ambiguous = ambiguous_vec))
                                                                              })
  # STEP 4: Writing matrices to disk ####
  
  names(rslt_btch) = cells_
  spliced_mat = unspliced_mat = ambiguous_mat = Matrix(data = 0, nrow = nrow(genes_tbl), ncol = length(cells_), sparse = T,
                                                       dimnames = list(gene = rownames(genes_tbl), cell = cells_))
  for(cell_ in cells_)
  {
    spliced_mat[,cell_]   = rslt_btch[[cell_]][,"spliced"]
    unspliced_mat[,cell_] = rslt_btch[[cell_]][,"unspliced"]
    ambiguous_mat[,cell_] = rslt_btch[[cell_]][,"ambiguous"]
  }
  
  null_ = writeMM(obj = spliced_mat,   file = paste0(btch_,'/','spliced.mtx'))      # sparse matrix for spliced counts
  null_ = writeMM(obj = unspliced_mat, file = paste0(btch_,'/','unspliced.mtx'))
  null_ = writeMM(obj = ambiguous_mat, file = paste0(btch_,'/','ambiguous.mtx'))
  cat(rownames(genes_tbl), sep = '\n', file = paste0(btch_,'/','genes.tsv'))        # genes
  cat(cells_, sep = '\n',              file = paste0(btch_,'/','barcodes.tsv'))     # cells
}

