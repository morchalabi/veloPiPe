library(readr)
library(parallel)

# STEP 1: Setting up options ####

inURL_ = 'ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz'      # URL+gtf file name
biotype_ = c('protein_coding')                                                                        # type of exons to extract
chrs_ = '^[1-9]+|^X|^Y|^MT'                                                                           # regular expression for extracting canonical and MT chromosomes
cpu_ = 3                                                                                              # number of CPU cores

# STEP 2: Reading in gtf file online ####
message('Reading in gtf file online')

mus_gtf = read_delim(file = inURL_, comment = '#', na = '.', delim = '\t',
                     col_names = c('seq','source','feature','start','end','score','strand','frame','attribute'), col_types = 'cccdddcic')

mus_gtf = mus_gtf[grepl(x = mus_gtf$seq, pattern = chrs_ , perl = T),c("seq","feature","start","end","attribute","strand")]     # extracting chromosomes/sequences
mus_gtf = mus_gtf[order(mus_gtf$seq, mus_gtf$start, mus_gtf$end),]                                                              # sorting gtf by seq, start and end

# STEP 3: Generating gene table ####
message('Generating gene table')

# not possible to filer out non-protein-coding genes here as some pseudo-genes have protein coding transcripts
genes_tbl = as.data.frame(mus_gtf[mus_gtf$feature %in% 'gene',])
genes_ = unlist(mclapply(X = genes_tbl$attribute, mc.cores = cpu_,  FUN =      # extracting all MGI genes
                                                                    function(att_)
                                                                    {
                                                                      annots_ = strsplit(x = att_, split = '; ')[[1]]
                                                                      annots_ = annots_[grepl(x = annots_, pattern = 'gene_name')]
                                                                      gene_ = strsplit(x = annots_[1], split = '"')[[1]][2]
                                                                      return(gene_)
                                                                    }))
inds_ = !duplicated(genes_, fromLast = T)     # keeps one copy of conflicting genes. Genes are sorted by start and then end; using fromLast = T make this command keep last observations having longer loci
genes_ = genes_[inds_]
genes_tbl = genes_tbl[inds_,c('seq','start','end')]
rownames(genes_tbl) = genes_

# STEP 4: Generating exon table ####
message('Generating exon table')

exons_tbl = mus_gtf[mus_gtf$feature %in% 'exon',]
gene_exon = do.call(mclapply(X = exons_tbl$attribute, mc.cores = cpu_,  FUN =      # extracting MGI protein-coding exons
                                                                        function(att_)
                                                                        {
                                                                          annots_ = strsplit(x = att_, split = '; ')[[1]]
                                                                          annots_ = annots_[grepl(x = annots_, pattern = 'gene_name|transcript_biotype|exon_id')]
                                                                          btyp_ = strsplit(x = annots_[2], split = '"')[[1]][2]
                                                                          if(btyp_ %in% biotype_)
                                                                          {
                                                                            gene_ = strsplit(x = annots_[1], split = '"')[[1]][2]
                                                                            exon_ = strsplit(x = annots_[3], split = '"')[[1]][2]
                                                                            return(c(gene = gene_, exon = exon_))
                                                                          }
                                                                          return(NA)
                                                                        }), what = rbind)

exons_tbl = cbind(gene_exon, exons_tbl[,c("seq","start","end", "strand")])
exons_tbl = exons_tbl[!is.na(exons_tbl$gene),]
exons_tbl = exons_tbl[!duplicated(exons_tbl[,c("start","end")]), ]      # one exon can be part of several transcripts
genes_tbl = genes_tbl[unique(exons_tbl$gene),]                          # some genes in genes_tbl are not protein coding

# STEP 5: Saving data ####

save(genes_tbl, exons_tbl, file = '../data/genes_exons_tables.RData')
