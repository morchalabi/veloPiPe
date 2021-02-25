library(parallel)

# STEP 1: Setting up options ####

inURI_ = '../data/'               # top directory of bam files
chrs_ = c(1:19,'X','Y','MT')      # chromosomes/sequences to extract from BAM files; current one is set for canonical chromosomes of mouse
MAPQ_ = 8                         # MAPQ = -10*log10(probability that the alignment is false)
FLAG_ = NULL                      # FLAG in samfile
cpu_ = 3                          # number of CPU cores

# STEP 2: Processing BAM files ####

bamfiles_ = list.files(path = inURI_, pattern = '*.bam\\b', full.names = T, recursive = T)
bamfiles_ = bamfiles_[!grepl(x = bamfiles_, pattern = '*.bai\\b')]

if(is.null(FLAG_))
{
  null_ = mclapply(X = bamfiles_, FUN =
                                  function(bamfile_)
                                  {
                                    out_bam = paste0(strsplit(bamfile_, split = '.bam')[[1]], '.processed.bam')
                                    system2(command = 'samtools', args = c('view', '-b', '-q', MAPQ_, bamfile_, chrs_, '-o', out_bam))                    # BAM file must be already position sorted
                                    system2(command = 'samtools', args = c('index', out_bam))
                                    return(NULL)
                                  }, mc.cores = cpu_)
}else
{
  null_ = mclapply(X = bamfiles_,  FUN =
                                   function(bamfile_)
                                   {
                                     out_bam = paste0(strsplit(bamfile_, split = '.bam')[[1]], '.processed.bam')
                                     system2(command = 'samtools', args = c('view', '-b', '-q', MAPQ_, '-f', FLAG_, bamfile_, chrs_, '-o', out_bam))      # BAM file must be already position sorted
                                     system2(command = 'samtools', args = c('index', out_bam))
                                     return(NULL)
                                   }, mc.cores = cpu_)
}
