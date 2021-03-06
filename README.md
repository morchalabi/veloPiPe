# veloPiPe
Universal count pipeline for RNA velocity estimation. VeloPiPe potentially can be applied to any scRNA-seq protocol including MARS-seq, SMART-seq, Dropt-seq, 10x Genomic, etc. as it only needs bam files containing final alignments. The pipline outputs spliced/unspliced matrices in Matrix Market and loom formats which are used as an input to RNA velocity estimation tools like [scVelo](https://scvelo.readthedocs.io/) and [velocyto](https://velocyto.org/).

![Anterior Primitive Streak cells](https://github.com/morchalabi/veloPiPe/blob/dev/doc/scVelo.jpg)
