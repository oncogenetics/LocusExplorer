*bedGraph* fields are:

`chrom` - The name of the chromosome (e.g. chr3, chrY).  
`chromStart` - The starting position of the feature in the chromosome.    
`chromEnd` - The ending position of the feature in the chromosome.  
`score` - A score, any number.    

See [BedGraph Track Format](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) for more details.

File is tab separated and has no header. This file will be used to create a bar chart. Score is the height, e.g.:

```
chr2	173292313	173371181	-100
chr2	173500000	173520000	1000
```

