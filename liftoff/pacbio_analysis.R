library(pafr)
library(seqinr)
library(stringr)
library(Biostrings)
library(intervals)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(plot.matrix)

#####
#LOAD PACBIO GENOME FASTA FILES AND REFERENCE GENOME
#####

fasta.files <- list.files("~/data/pacbio RR genomes/input/fasta files/")[seq(3, 33, by = 2)]
fastas <- list(16)
for (x in 1:16) fastas[[x]] <- read.fasta(paste0("~/data/pacbio RR genomes/input/fasta files/", fasta.files[x], collapse = "")) 
names(fastas) <- sapply(24:39, function (x) paste0("MSY", x, collapse = ""))

ref.genome <- read.fasta("~/data/pacbio RR genomes/input/SacCer3.fna")
chrom.lengths <- sapply(ref.genome, length) 

chrom.name.table <- cbind.data.frame(stringsAsFactors = F,
                                     strain = rep(sapply(24:39, function (x) paste0("MSY", x, collapse = "")), each = 17),
                                     contig = unlist(lapply(fastas, getName)),
                                     chrom = unlist(lapply(1:16, function (x) sapply(getAnnot(fastas[[x]]), function (y) strsplit(y, ", ")[[1]][2]))))

chrom.name.table$chrom[which(chrom.name.table$chrom == "mitochondrion")] <- "chrmt"

chrom.name.list <- split(chrom.name.table, chrom.name.table$strain)

######
#FUNCTIONS TO GET SPECIFIC DNA SEQUENCES
######

dna.sequence.from.file <- function (strain, chrom, start, end, rev.comp = F, directory = "~/data/pacbio RR genomes/input/fasta files/") {
  files <- list.files(directory)
  fasta.file <- paste0(directory, files[intersect(grep(strain, files, ignore.case = T),grep("fasta$", files))], collapse = "")
  fasta.seq <- read.fasta(fasta.file)
  dna.seq <- fasta.seq[[chrom]][start:end]
  if (rev.comp) dna.seq <- rev(comp(dna.seq))
  paste0(dna.seq, collapse = "")
}

dna.sequence <- function (strain.fasta, chrom, start, end, rev.comp=F) {
  if (start < 1 || end > length(strain.fasta[[chrom]])) "" else {
    dna.seq <- strain.fasta[[chrom]][start:end]
    if (rev.comp) dna.seq <- rev(comp(dna.seq))
    paste0(dna.seq, collapse = "")
  }
}

dna.sequence.from.gff <- function (strain,
                                   gff.line,
                                   padding = 0, padding.l = padding, padding.r = padding,
                                   gff.file = gff.files[[grep(strain, names(gff.files), ignore.case = T)]],
                                   fasta.file = fastas[[grep(strain, names(fastas), ignore.case = T)]]) {
  if (gff.file$V7[gff.line] == "+") rc <- F else rc <- T
  dna.sequence(fasta.file, gff.file$V1[gff.line], as.numeric(gff.file$V4[gff.line]) - padding.l, as.numeric(gff.file$V5[gff.line]) + padding.r, rev.comp = rc)
}

#####
#LOAD QTL info determined by Josh (eLife 2019)
#####

#LOD scores and other QTL info from Josh, sent to me prior to publication of his eLife 2019 paper. I used these to tell which 
#QTLs were statistically significant, but not their LOD scores, which I calculated from the new plots I made (I'm not sure why 
#precisely they differ slightly - his LOD scores are higher than what I calculate. I know he did some filtering of poorly-behaving segregants).
qtl.table <- read.csv("~/data/pacbio RR genomes/input/QTL_data/round_robin_QTLs.csv", stringsAsFactors = F, header = T)
crosses <- c("375", "A", "376", "B", "377", "393", "381", "3008", "2999", "3000", "3001", "3049", "3003", "3004", "3043", "3028")
for (x in 1:15) qtl.table$cross[which(qtl.table$cross == crosses[x])] <- paste0("MSY",x+23, "_","MSY", x+24, collapse = "")
qtl.table$cross[which(qtl.table$cross == crosses[16])] <- "MSY39_MSY24"

#The QTL data from the published dataset (above is QTL data I obtained from Josh directly pre-publication). 
#Available here: https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDkyMTIvZWxpZmUtNDkyMTItZmlnMy1kYXRhMS12Mi54bHM-/elife-49212-fig3-data1-v2.xls?_hash=2eBqo5QmCeNhn%2B%2FQIUQPcgFMrsc68eUCmEecg%2Fp3O0I%3D
qtl.table.published <- read.csv("~/data/pacbio RR genomes/input/QTL_data/elife-49212-fig3-data1-v2.csv", stringsAsFactors = F, header = T)
nrow(qtl.table.published) ; nrow(qtl.table.published)/(16*38) #7751 QTLs mapped in the published dataset; 12.7 per trait per cross.
#Median size of 31.7 kb in the published dataset (1.5 LOD drop)
median(sapply(qtl.table.published$X1.5.LOD.drop.CI..right, function (x) as.numeric(strsplit(x, "_")[[1]][2])) - sapply(qtl.table.published$X1.5.LOD.drop.CI..left, function (x) as.numeric(strsplit(x, "_")[[1]][2])))

phenos <- read.delim("~/data/pacbio RR genomes/input/QTL_data/phenotypes.tsv") #I got this from https://github.com/joshsbloom/yeast-16-parents

##########
#GENE FAMILY ANNOTATION
##########

#I got the gene family dataset of Genolevures from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099480, file S1
genefam <- scan("~/data/pacbio RR genomes/input/pone.0099480.s001.txt", what = "character")

genefamlist <- list()
#There are an equal number of "cerevisiae" and "glabrata", and "glabrata" always comes after "cerevisiae".
#But sometimes there are gene families missing from cerevisiae entirely.
cerind <- grep("cerevisiae", genefam)
glaind <- grep("glabrata", genefam) 
cerind.nonzero <- cerind[-which((cerind + 1) %in% glaind)]
glaind.nonzero <- glaind[-which((cerind + 1) %in% glaind)]
for (x in 1:length(cerind.nonzero)) genefam[(cerind.nonzero[x] + 1):(glaind.nonzero[x] - 1)] -> genefamlist[[x]]

names(genefamlist) <- genefam[cerind.nonzero - 1]

#Getting just the systematic yeast gene name
genefamlist2 <- lapply(genefamlist, function (x) sapply(x, function (y) strsplit(y, "\\(")[[1]][2]))
genefamlist2 <- lapply(genefamlist2, function (x) sapply(x, function (y) strsplit(y, "\\)")[[1]][1]))

#Some genes are missing their systematic names, for some reason. Adding them back in:
nafams <- c("GL3C0059",
            "GL3C0702",
            "GL3C2818",
            "GL3C3871",
            "GL3C3871",
            "GL3C3108",
            "GL3C3108",
            "GL3C2818",
            "GL3C4502",
            "GL3C3785",
            "GL3R2024",
            "GL3C3841",
            "GL3C3841",
            "GL3C4314",
            "GL3R1741")
nagenes <- c(2, 1, 1, 1, 2, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1)
extrasystematicnames <- c("YBL112C",
                          "YBR016W",
                          "YCL048W-A",
                          "YDL184C",
                          "YDL134C",
                          "YDR034W-B",
                          "YDR210W",
                          "YDR524C-B",
                          "YER188W",
                          "YKL096W-A",
                          "YLR262C-A",
                          "YMR251W-A",
                          "YOL052C-A",
                          "YOR302W",
                          "YPR169W-A")

for (x in 1:15) genefamlist2[[nafams[x]]][nagenes[x]] <- extrasystematicnames[x]

genefamtable <- cbind.data.frame(genes = unlist(genefamlist2), fams = rep(names(genefamlist2), times=lengths(genefamlist2)))

save(genefamtable, file = "~/data/pacbio RR genomes/output/genefamtable.R")
save(genefamlist2, file = "~/data/pacbio RR genomes/output/genefamlist.R")
load("~/data/pacbio RR genomes/output/genefamtable.R")
load("~/data/pacbio RR genomes/output/genefamlist.R")

#######
#GFF FILES
#######

#A function to add gene family identifiers (from above) to a gff file object.
add.gene.fams <- function (gff, genefamilytable = genefamtable) {
  gff <- cbind.data.frame(gff, ID = sapply(gff$V9, function (x) strsplit(strsplit(x, "Name=")[[1]][2], ";")[[1]][1]), stringsAsFactors = F)
  genefam.assignments <- sapply(gff[,10], function (x) unique(genefamtable$fams[which(genefamtable$genes == as.character(x))]))
  genefam.assignments <- sapply(genefam.assignments, function (x) if (length(x)) as.character(x) else NA)
  gff <- cbind.data.frame(gff, GENE.FAM = genefam.assignments, stringsAsFactors = F)
  gff
}

#Read in the gff files made earlier by Liftoff on biowulf, by lifting over gene annotations from the reference genome onto the pacbio genomes.
gff.files <- list()
for (x in 1:16) gff.files[[x]] <- read.delim(header = F, stringsAsFactors = F, paste0("~/data/pacbio RR genomes/input/gff_files/MSY", x + 23, ".gff", collapse = ""), quote = "", fill = F)
names(gff.files) <- sapply(24:39, function (x) paste0("msy", x, "gff", collapse = ""))

gff.files.fams <- lapply(gff.files, add.gene.fams)

genefamcount <- matrix(nrow = length(unique(genefamtable$fams)), ncol = 16)
colnames(genefamcount) <- sapply(24:39, function (x) paste0("MSY", x, collapse = ""))
rownames(genefamcount) <- unique(genefamtable$fams)
for (j in 1:16){
  for (i in 1:length(unique(genefamtable$fams))) {
    genefamcount[i,j] <- length(which(gff.files.fams[[j]]$GENE.FAM == rownames(genefamcount)[i]))
  }
}

save(genefamcount, file = "~/data/pacbio RR genomes/output/genefamcount.R")
load("~/data/pacbio RR genomes/output/genefamcount.R")

for (x in 1:16) {
  chrom.temp <- sapply(gff.files[[x]]$V1, function (y) chrom.name.list[[x]]$chrom[which(chrom.name.list[[x]]$contig == y)])
  gff.files[[x]] <- cbind.data.frame(gff.files[[x]],
                                     CHROM = chrom.temp)
  gff.files.fams[[x]] <- cbind.data.frame(gff.files.fams[[x]],
                                          CHROM = chrom.temp)
}
for (x in 1:length(gff.files)) rownames(gff.files[[x]]) -> rownames(gff.files.fams[[x]])

save(gff.files, file = "~/data/pacbio RR genomes/output/gff.files.R")
save(gff.files.fams, file = "~/data/pacbio RR genomes/output/gff.files.fam.R")
load("~/data/pacbio RR genomes/output/gff.files.fam.R")
load("~/data/pacbio RR genomes/output/gff.files.R")

###########
#COMBINING DE NOVO ORF FINDING INTO THE LIFTOFF-GENERATED GFFs
###########

#Using ORF-finding data done by Michael
pacbio.blast <- read.csv("~/data/pacbio RR genomes/input/michael_ORF_finding/all-hits_s288c-orfs.csv", stringsAsFactors = F)

#Note that Michael's ORFs have a slightly off indexing: the A of the ATG starts from the base after the "start" coordinate, as well as the final base of the final codon, whereas for other genes we have it
#with the A starting from the start coordinate, and likewise for the final codon. I will leave it as is. If we change it later, remember to keep the start and end coordinates within the chromosomes - I think
#some of Michael's ORFs go right to the chromosome end.

percent_id_whole <- sapply(1:nrow(pacbio.blast),function (x) pacbio.blast$length[x]*3/(abs(pacbio.blast$start[x] - pacbio.blast$end[x]) - 2) * pacbio.blast$percent_id[x]/100)
pacbio.blast <- cbind.data.frame(pacbio.blast, percent_id_whole)

pacbio.blast$name <- paste(pacbio.blast$name, pacbio.blast$chr, sep = "_")

pacbio.blast.top.hits <- do.call("rbind", lapply(split(pacbio.blast, pacbio.blast$name), function (x) x[order(x$percent_id_whole, decreasing = T)[1],]))

pacbio.no.blast <- read.csv("~/data/pacbio RR genomes/input/michael_ORF_finding/unique-from-s288c.csv", stringsAsFactors = F)[-268,] #Row 268 is labelled as coming from "S288C", and I don't know if that means MSY25 or the reference genome. Doesn't look like a real gene, in any case

pacbio.no.blast$name <- paste(pacbio.no.blast$name, pacbio.no.blast$chr, sep = "_")

pacbio.blast.top.hits$chr[which(pacbio.blast.top.hits$chr == "mitochondrion")] <- "chrmt"
pacbio.no.blast$chr[which(pacbio.no.blast$chr == "mitochondrion")] <- "chrmt"

pacbio.orfs.gff <- cbind.data.frame(V2 = rep("ORF-mining", nrow(pacbio.blast.top.hits)),
                                    V3 = rep("ORF", nrow(pacbio.blast.top.hits)),
                                    V4 = apply(pacbio.blast.top.hits, 1, function (x) min(x[3:4])),
                                    V5 = apply(pacbio.blast.top.hits, 1, function (x) max(x[3:4])),
                                    V6 = rep(".", nrow(pacbio.blast.top.hits)),
                                    V7 = apply(pacbio.blast.top.hits, 1, function (x) if (x[3] < x[4]) "+" else "-"),
                                    V8 = rep(".", nrow(pacbio.blast.top.hits)),
                                    V9 = apply(pacbio.blast.top.hits, 1, function (x) paste0("ID=", x[1], ";Note=BLAST hit with ", as.numeric(x[7]), " protein percent ID to ", x[6], " across ", round(300*as.numeric(x[8])/(abs(as.numeric(x[3]) - as.numeric(x[4])) - 2), digits = 3), " percent of the ORF length, with e-value ", as.numeric(x[10]))),
                                    CHROM = pacbio.blast.top.hits$chr)

pacbio.orfs.gff <- rbind.data.frame(pacbio.orfs.gff, 
                                    cbind.data.frame(V2 = rep("ORF-mining", nrow(pacbio.no.blast)),
                                                     V3 = rep("ORF", nrow(pacbio.no.blast)),
                                                     V4 = apply(pacbio.no.blast, 1, function (x) min(x[3:4])),
                                                     V5 = apply(pacbio.no.blast, 1, function (x) max(x[3:4])),
                                                     V6 = rep(".", nrow(pacbio.no.blast)),
                                                     V7 = apply(pacbio.no.blast, 1, function (x) if (x[3] < x[4]) "+" else "-"),
                                                     V8 = rep(".", nrow(pacbio.no.blast)),
                                                     V9 = apply(pacbio.no.blast, 1, function (x) paste0("ID=", x[1], ";Note=No BLAST hits to s288c, with sequence ", x[6])),
                                                     CHROM = pacbio.no.blast$chr))

#Use chrom.name.list to add the V1 column  
chrom.name.table.short <- chrom.name.table
chrom.name.table.short$chrom <- sapply(chrom.name.table.short$chrom, function (x) strsplit(x, " |-")[[1]][1]) #A quick check indicated Michael only used the starts of the longer descriptive chromosome names.
chrom.name.list.short <- split(chrom.name.table.short, chrom.name.table.short$strain)
chrom.name.list.short$MSY31$chrom[c(4, 13)] <- c("chrV", "chrXIV") #Manually correcting a annotation mistake caused by a translocation

V1 <- apply(pacbio.blast.top.hits, 1, function (x) chrom.name.list.short[[substr(x[1], 1, 5)]]$contig[which(chrom.name.list.short[[substr(x[1], 1, 5)]]$chrom == x[2])])
V1 <- c(V1, 
        apply(pacbio.no.blast, 1, function (x) chrom.name.list.short[[substr(x[1], 1, 5)]]$contig[which(chrom.name.list.short[[substr(x[1], 1, 5)]]$chrom == x[2])]))

pacbio.orfs.gff <- cbind.data.frame(V1, pacbio.orfs.gff)

pacbio.orfs.gff.by.strain <- split(pacbio.orfs.gff, sapply(pacbio.orfs.gff$V9, function (x) substr(x, 4, 8)))

#Avoid ORFs that overlap other extra elements than genes (also included Y' elements, X elements, transposons, centromeres, tRNAs, and ARSes)
extra.gff.files <- list()
for (x in 1:16) extra.gff.files[[x]] <- read.delim(header = F, stringsAsFactors = F, paste0("~/data/pacbio RR genomes/input/gff_files/221019 - extra/MSY", x + 23, "extra.gff", collapse = ""), quote = "", fill = F)
names(extra.gff.files) <- sapply(24:39, function (x) paste0("msy", x, "gff", collapse = ""))

#Pulling out ORFs that are in gaps between Liftoff genes and adding them to the gff files
extra.orfs <- list()

for (x in 1:16) {
  extra.orfs[[x]] <- do.call("rbind", lapply(unique(extra.gff.files[[x]]$V1), function (y) pacbio.orfs.gff.by.strain[[x]][which(pacbio.orfs.gff.by.strain[[x]]$V1 == y),][which(lengths(interval_overlap(Intervals(sapply(pacbio.orfs.gff.by.strain[[x]][which(pacbio.orfs.gff.by.strain[[x]]$V1 == y), 4:5], as.numeric)), Intervals(extra.gff.files[[x]][which(extra.gff.files[[x]]$V1 == y), c("V4", "V5")]))) == 0), ]))
}

extra.orfs.over.100codons <- lapply(extra.orfs, function (x) x[which((as.numeric(x$V5) - as.numeric(x$V4)) > 300),])

combined.gffs <- list()

for (x in 1:16) {
  combined.gffs[[x]] <- rbind.data.frame(gff.files[[x]], extra.orfs.over.100codons[[x]])
}

names(combined.gffs) <- names(gff.files)

for (x in 1:16) {
  combined.gffs[[x]] <- combined.gffs[[x]][order(combined.gffs[[x]]$V1, as.numeric(combined.gffs[[x]]$V4)),]
}

#It turns out a lot of the "ORFs" don't have a start codon. The theory is that they're translating from alternate start codons, but I don't believe that's likely. Removing them below:
for (i in 1:16) {
  orfaas <- translate(DNAStringSet(sapply(which(combined.gffs[[i]]$V3 == "ORF"), function (x) dna.sequence.from.gff(paste0("MSY", 23 + i), gff.line = x, gff.file = combined.gffs[[i]], padding.l = -1, padding.r = 1))), no.init.codon = T)
  orfaas
  orfswoatg <- which(sapply(orfaas, function (x) !grepl("M", substr(x, 1, nchar(x) - 101))))
  combined.gffs[[i]] <- combined.gffs[[i]][-(which(combined.gffs[[i]]$V3 == "ORF")[orfswoatg]),]
}

#Writing the combined GFF files
strains <- c("M22", "BY","RM","YPS163","YJM145","CLIB413","YJM978","YJM454","YPS1009","I14","Y10","PW5","273614","YJM981","CBS2888", "CLIB219")
for (x in 1:16) {
  write.table(combined.gffs[[x]], file = paste0("~/data/pacbio RR genomes/output/gffs/", strains[x], " (MSY", x + 23, ").gff"), quote = FALSE, sep = "\t", row.names = F, col.names = F)
}

combined.gffs <- list()
for (i in 1:16) combined.gffs[[i]] <- read.delim(paste0("~/data/pacbio RR genomes/output/gffs/", grep(paste0("MSY", i + 23), list.files("~/data/pacbio RR genomes/output/gffs/"), value = T)), header = F, stringsAsFactors = F, quote = "", fill = F)
names(combined.gffs) <- names(gff.files)

#######
#DETERMINING WHETHER THE END OF THE ASSEMBLED CHROMOSOMES HAVE THE COMPONENTS THAT ENDS OF CHROMOSOMES ARE EXPECTED TO HAVE
#######

#Table of chromosome ends to contigs
chrom.end.content <- data.frame(matrix(ncol = 8, nrow = 512))
colnames(chrom.end.content) <- c("strain", "contig", "contig.end", "ref.chrom.equivalent", "ref.chrom.end", "telomeric.repeats", "x_element", "y_element")
chrom.end.content$strain <- rep(paste0("MSY", 24:39), each = 32)
chrom.end.content$ref.chrom.equivalent <- rep(rep(1:16, each = 2), 16)
chrom.end.content$ref.chrom.end <- rep(c("L", "R"), 256)

#Most of the contigs were associated to a single chromosome by Morgan Park
for (i in 1:16) {
  for (j in 1:17) {
    if (nchar(chrom.name.list[[i]]$chrom[j]) < 20 && chrom.name.list[[i]]$chrom[j] != "chrmt") {
      chrom.end.content$contig[intersect(which(chrom.end.content$strain == paste0("MSY", i + 23)),
                                         which(chrom.end.content$ref.chrom.equivalent == as.numeric(as.roman(substr(chrom.name.list[[i]]$chrom[j], 4, nchar(chrom.name.list[[i]]$chrom[j]))))))] <- chrom.name.list[[i]]$contig[j]
    }
  }
}

for (i in 1:512) {
  if (!is.na(chrom.end.content$contig[i])) chrom.end.content$contig.end[i] <- c("R", "L")[1 + (i %% 2)]
}

#A few contigs were formed by reciprocal translocations between chromosomes. Here I manually annotated which chromosome end is associated to which contig end for those cases. To do that I used alignments to the reference genome made by minimap2
chrom.end.content[c(15, 16, 31, 32), c("contig", "contig.end")] <- cbind.data.frame(c("tig00000012rc", "tig00000006", "tig00000006", "tig00000012rc"), c("L", "R", "L", "R"))
chrom.end.content[c(233, 234, 251, 252), c("contig", "contig.end")] <- cbind.data.frame(c("tig00000009rc", "tig00000030rc", "tig00000009rc", "tig00000030rc"), c("R", "L", "L", "R"))
chrom.end.content[c(431 ,432, 447, 448), c("contig", "contig.end")] <- cbind.data.frame(c("tig00000134rc", "tig00000120", "tig00000120", "tig00000134rc"), c("L", "R", "L", "R"))
chrom.end.content[c(467, 468, 477, 478), c("contig", "contig.end")] <- cbind.data.frame(c("tig00000004", "tig00000020rcj44", "tig00000020rcj44", "tig00000004"), c("L", "R", "L", "R"))

#Telomeric repeats
for (i in seq(1, 512, by = 2)) {
  chrom.end.content$telomeric.repeats[i] <- str_count(paste0(head(fastas[[chrom.end.content$strain[i]]][[chrom.end.content$contig[i]]], 100), collapse = ""), "ac")
}

for (i in seq(2, 512, by = 2)) {
  chrom.end.content$telomeric.repeats[i] <- str_count(paste0(tail(fastas[[chrom.end.content$strain[i]]][[chrom.end.content$contig[i]]], 100), collapse = ""), "tg")
}

#X elements
for (i in 24:39) {
  x_element_gff <- extra.gff.files[[grep(i, names(extra.gff.files), value = T)]][grep("X_element", extra.gff.files[[grep(i, names(extra.gff.files), value = T)]]$V9), 1:8]
  for (j in 1:32) {
    contig <- chrom.end.content$contig[grep(i, chrom.end.content$strain)][j]
    arm <- chrom.end.content$contig.end[grep(i, chrom.end.content$strain)][j]
    if (contig %in% x_element_gff$V1) {
      if (arm == "L") {
        if (min(x_element_gff$V4[which(x_element_gff$V1 == contig)]) < 30000) chrom.end.content$x_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                              which(chrom.end.content$contig == contig)),
                                                                                                                    which(chrom.end.content$contig.end == "L"))] <- T
      } else {
        if (max(x_element_gff$V4[which(x_element_gff$V1 == contig)]) > 200000) chrom.end.content$x_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                               which(chrom.end.content$contig == contig)),
                                                                                                                     which(chrom.end.content$contig.end == "R"))] <- T
      }
    } 
  }
}

#Y' elements
for (i in 24:39) {
  y_element_gff <- extra.gff.files[[grep(i, names(extra.gff.files), value = T)]][grep("Y_prime_element", extra.gff.files[[grep(i, names(extra.gff.files), value = T)]]$V9), 1:8] #These are annotations of the y' element itself.
  alt_y_element_gff <- gff.files.fams[[grep(i, names(extra.gff.files), value = T)]][which(gff.files.fams[[grep(i, names(extra.gff.files), value = T)]]$GENE.FAM == "GL3C0059"), 1:8] #There are also annotations of the y' helicase as a gene
  for (j in 1:32) {
    contig <- chrom.end.content$contig[grep(i, chrom.end.content$strain)][j]
    arm <- chrom.end.content$contig.end[grep(i, chrom.end.content$strain)][j]
    if (contig %in% y_element_gff$V1) {
      if (arm == "L") {
        if (min(y_element_gff$V4[which(y_element_gff$V1 == contig)]) < 30000) chrom.end.content$y_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                              which(chrom.end.content$contig == contig)),
                                                                                                                    which(chrom.end.content$contig.end == "L"))] <- T
      } else {
        if (max(y_element_gff$V4[which(y_element_gff$V1 == contig)]) > 200000) chrom.end.content$y_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                               which(chrom.end.content$contig == contig)),
                                                                                                                     which(chrom.end.content$contig.end == "R"))] <- T
      }
    }
    if (contig %in% alt_y_element_gff$V1) {
      if (arm == "L") {
        if (min(alt_y_element_gff$V4[which(alt_y_element_gff$V1 == contig)]) < 30000) chrom.end.content$y_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                                      which(chrom.end.content$contig == contig)),
                                                                                                                            which(chrom.end.content$contig.end == "L"))] <- T
      } else {
        if (max(alt_y_element_gff$V4[which(alt_y_element_gff$V1 == contig)]) > 200000) chrom.end.content$y_element[intersect(intersect(grep(i, chrom.end.content$strain),
                                                                                                                                       which(chrom.end.content$contig == contig)),
                                                                                                                             which(chrom.end.content$contig.end == "R"))] <- T
      }
    }  
  }
}


plot(xlim = c(0, 130), ylim = c(80, 40), type = "n", 1, axes = F, ylab = "", xlab = "chromosome", asp = 1)

axis(side = 1, at = seq(4, 130, by = 8), labels = as.roman(1:16), tick = F, cex.axis = 0.7)
axis(side = 2, at = seq(48, 78, by = 2), labels = strains, tick = F, las = 2, cex.axis = .7)

for (i in 1:512) {
  ypos <- as.numeric(substr(chrom.end.content$strain[i], 4, 5)) * 2
  lr <- ifelse(chrom.end.content$ref.chrom.end[i] == "L", -1, 1) #whether we're looking at the left or right end of the chromosome
  xpos <- chrom.end.content$ref.chrom.equivalent[i] * 8 - 4.5 #Add 0.5 space more than needed so there is a little space between the L and R triplet.
  lines(c(xpos + 0.25, xpos + 0.75), c(ypos, ypos), col = "#606060")
  polygon(c(xpos + lr * 2.75, xpos + lr * 2.75, xpos + lr * 2.75 + 1, xpos + lr * 2.75 + 1), c(ypos - 0.45, ypos + 0.45, ypos + 0.45, ypos - 0.45), col = ifelse(chrom.end.content$telomeric.repeats[i] > 20, "#333333", "white"), border = "#333333") #left telomere repeats
  polygon(c(xpos + lr * 1.75, xpos + lr * 1.75, xpos + lr * 1.75 + 1, xpos + lr * 1.75 + 1), c(ypos - 0.45, ypos + 0.45, ypos + 0.45, ypos - 0.45), col = ifelse(chrom.end.content$y_element[i] == T, "#666666", "white"), border = "#333333") #y' elements. For color, could use magenta
  polygon(c(xpos + lr * 0.75, xpos + lr * 0.75, xpos + lr * 0.75 + 1, xpos + lr * 0.75 + 1), c(ypos - 0.45, ypos + 0.45, ypos + 0.45, ypos - 0.45), col = ifelse(chrom.end.content$x_element[i] == T, "#999999", "white"), border = "#333333") #x' elements. For color, could use gold.
}

#########
#KNOWN MAJOR CHROMOSOMAL REARRANGEMENTS 
#########

#3 of our strains (M22, YJM981, and YJM454) have previously been found to have translocations (PMID 24814147). Below, I checked whether the genes flanking the
#known translocations are indeed found on the contigs we'd expect based on the translocation.

#Compare the chromosomes these genes are found on to fig. 4 from PMID 24814147
#Strain M22 has a chrXVI/VIII translocation
gff.files$msy24gff[grep("ECM34", gff.files$msy24gff$V9), -9] #chrVIII in ref genome; contig tig00000006 here
gff.files$msy24gff[grep("NOG1", gff.files$msy24gff$V9), -9] #chrXVI in ref genome; contig tig00000006 here
gff.files$msy24gff[grep("SSU1", gff.files$msy24gff$V9), -9] #chrXVI in ref genome; contig tig00000012rc here
gff.files$msy24gff[grep("YHL044W", gff.files$msy24gff$V9), -9] #chrVIII in ref genome; contig tig00000012rc here
gff.files$msy24gff[grep("GLR1", gff.files$msy24gff$V9), -9] #chrXVI in ref genome; contig tig00000012rc here

#Strain YJM981 has the same chrXVI/VIII translocation
gff.files$msy37gff[grep("ECM34", gff.files$msy37gff$V9), -9] #chrVIII in ref genome; contig tig00000120  here
gff.files$msy37gff[grep("NOG1", gff.files$msy37gff$V9), -9] #chrXVI in ref genome; contig tig00000120  here
gff.files$msy37gff[grep("SSU1", gff.files$msy37gff$V9), -9] #chrXVI in ref genome; contig tig00000134rc here
gff.files$msy37gff[grep("YHL044W", gff.files$msy37gff$V9), -9] #chrVIII in ref genome; contig tig00000134rc here
gff.files$msy37gff[grep("GLR1", gff.files$msy37gff$V9), -9] #chrXVI in ref genome; contig tig00000134rc here

#Strain YJM454 has a chrV/XIV translocation
gff.files$msy31gff[grep("YER132C", gff.files$msy31gff$V9), -9] #PMD1. chrV in ref genome; contig tig00000009rc here
gff.files$msy31gff[grep("PHO23", gff.files$msy31gff$V9), -9] #chrXIV in ref genome; contig tig00000009rc here
gff.files$msy31gff[grep("YNL098C", gff.files$msy31gff$V9), -9] #RAS2. chrXIV in ref genome; contig tig00000009rc here
gff.files$msy31gff[grep("YNL095C", gff.files$msy31gff$V9), -9] #chrXIV in ref genome; contig tig00000030rc here
gff.files$msy31gff[grep("GLC7", gff.files$msy31gff$V9), -9] #chrV in ref genome; contig tig00000030rc here

#Linkage analysis of a fourth translocation in YJM981 was done separately by Cory. 

#####
#COMPARING GENE CONTENT OF CHROMOSOMES
#####

#For each strain, make a data.frame of their gff genes with the following info: ID, orf_classification, sequence_ID, extra_copy_number

gene.table <- list()

for (i in 1:16) {
  gene.table[[i]] <- data.frame(
    ID = unname(sapply(gff.files[[i]]$V9[which(gff.files[[i]]$V3 == "gene")], function (x) strsplit(x, "ID=|;")[[1]][2])),
    orf_classification = unname(sapply(gff.files[[i]]$V9[which(gff.files[[i]]$V3 == "gene")], function (x) strsplit(strsplit(x, "orf_classification=")[[1]][2], ";")[[1]][1])),
    sequence_ID = unname(sapply(gff.files[[i]]$V9[which(gff.files[[i]]$V3 == "gene")], function (x) as.numeric(strsplit(strsplit(x, "sequence_ID=")[[1]][2], ";")[[1]][1]))),
    coverage = unname(sapply(gff.files[[i]]$V9[which(gff.files[[i]]$V3 == "gene")], function (x) as.numeric(strsplit(strsplit(x, "coverage=")[[1]][2], ";")[[1]][1]))),
    extra_copy_number = unname(sapply(gff.files[[i]]$V9[which(gff.files[[i]]$V3 == "gene")], function (x) as.numeric(strsplit(strsplit(x, "extra_copy_number=")[[1]][2], ";")[[1]][1])))
  )
}

for (i in 1:16) {
  gene.table[[i]] <- cbind.data.frame(gene.table[[i]], sequence_id_in_covered_region = gene.table[[i]]$sequence_ID/gene.table[[i]]$coverage)
}

gene.classifications <- list()

for (i in 1:16) {
  gene.classifications[[i]] <- apply(gene.table[[i]], 1, function (x) ifelse(x[2] == "Dubious", "dubious", 
                                                                             ifelse(as.numeric(x[5]) == 0 && as.numeric(x[6]) < 0.95, "diverged", 
                                                                                    ifelse(as.numeric(x[5]) == 0 && as.numeric(x[6]) >= 0.95, "retained", 
                                                                                           ifelse(as.numeric(x[5]) > 0 && as.numeric(x[6]) >= 0.95, "additional.copy", 
                                                                                                  ifelse(as.numeric(x[5]) > 0 && as.numeric(x[6]) < 0.95, "additional.homolog", NA))))))
}

ref.gff <- read.delim(header = F, stringsAsFactors = F, "~/data/pacbio RR genomes/input/sacCer3_sgd_short.gff", quote = "", fill = F)
ref.gff <- as.data.frame(matrix(ref.gff[22:160365,1], ncol = 9, byrow = T), stringsAsFactors = F)
ref.gff.gene.IDs <- sapply(ref.gff$V9[which(ref.gff$V3 == "gene")], function (x) strsplit(x, "=|;")[[1]][2])
ref.gff.gene.IDs <- unname(ref.gff.gene.IDs)

#lost genes (missing from table), diverged genes (cn = 1, homology 30-95%), retained genes (cn = 1, homology > 95%), additional copies (cn > 1, homology > 95%), additional diverged homologs (cn > 1, homology 30-95%), 
#additional ORFs (orfFinder).

gene.classification.table <- rbind(sapply(gene.classifications, table), 
                                   ORFs = sapply(combined.gffs, function (x) length(which(x$V3 == "ORF"))), 
                                   absent = sapply(1:16, function (x) length(which(is.na(match(ref.gff.gene.IDs[-grep("orf_classification=Dubious", ref.gff$V9[which(ref.gff$V3 == "gene")])], gene.table[[x]]$ID))))))

#Make a barplot. I think make it first without the negative values but extend the ylim to c(-150, 6500) and then add in the negative values with add = T.

barplot(gene.classification.table[c("retained", "diverged", "additional.copy", "additional.homolog", "ORFs"),], ylim = c(-max(gene.classification.table[7,]), max(apply(gene.classification.table[-7,], 2, sum))), col = c("#CCCCCC", "gold", "forestgreen", "cyan", "#EAB3F3"), border = NA, ylab = "Count", axes = F, xaxt = 'n')
barplot(add = T, -gene.classification.table[7,], col = "darkorange", border = NA, axes = F, xaxt = 'n')
axis(2, col = NA, col.ticks = "#606060")
axis(1, las = 2, at = seq(.7, .7 + 1.2 * 15, by = 1.2), strains, tick = F)
legend(4, 3000, c("Retained reference genes", "Genes diverged from reference", "Additional reference gene copies", "Additional diverged gene copies", "Additional ORFs", "Absent reference genes"),
       fill = c("#CCCCCC", "gold", "forestgreen", "cyan", "#EAB3F3", "darkorange"),
       border = NA, box.lwd = 0)

#########
#MALTOSE AND PARAQUAT LOD PLOTS
#########

lods <- function (phenos, geno.table) {
  as.vector(-length(phenos) * log(1 - cor(phenos, geno.table[,-1], use = 'pairwise.complete.obs')^2/(2*log(10))))
}

lod.plot <- function (lod.vector, geno.markers, highlight.markers = NA, title = NA, col1 = "cornflowerblue", col2 = "blue", lwd = 2, chrom = "genome", coords = NULL, ymax = 0, add = F, y.offset = 0) {
  drop.values <- which(is.na(lod.vector))
  chroms <- sapply(strsplit(geno.markers, "_"), function (x) x[1]) #This used to be [-c(1, drop.values)]; I'm not sure why the first element was dropped? Same for "positions" below.
  positions <- as.numeric(sapply(strsplit(geno.markers, "_"), function (x) x[2]))
  if (length(drop.values) > 0) { lod.vector <- lod.vector[-drop.values] ; chroms <- chroms[-drop.values] ; positions <- positions[-drop.values] }
  if (chrom == "genome") {
    if (add == F) {
      plot(1, type = "n", xlim = c(0, sum(chrom.lengths[-17]) + 80000*15), ylim = c(0, max(c(lod.vector, ymax))), axes = F, xlab = "chrom", ylab = "LOD", main = title)
      axis(2, col = NA, col.ticks = 1)
      for (x in 1:16) {
        chrom.lims <- c(seq(0,80000*16, by= 80000) + cumsum(c(0,chrom.lengths))[1:17])[x:(x + 1)] + c(0, -80000)
        axis(1, at = seq(chrom.lims[1], chrom.lims[2], length.out = 3), 
             labels = as.roman(c(NA, x, NA)), 
             col = "#777777", lwd.ticks = 0, lwd = 2, lend = 1)
      }
    }
    sapply(1:16, function (x) lines(positions[which(chroms == names(chrom.lengths)[x])] + 80000*(x-1) + cumsum(c(0,chrom.lengths))[x],
                                    lod.vector[which(chroms == names(chrom.lengths)[x])] + y.offset,
                                    col = c(col2, col1)[1 + (x %% 2)],
                                    lwd = lwd, lend = 1))
    if (!is.na(highlight.markers[1])) {
      highlight.chroms <- sapply(strsplit(highlight.markers, "_"), function (x) x[1])
      highlight.chrom.nums <- sapply(highlight.chroms, function (x) which(names(chrom.lengths) == x))
      highlight.positions <- as.numeric(sapply(strsplit(highlight.markers, "_"), function (x) x[2]))
      for (x in 1:length(highlight.markers)) abline(col = "red", lty = 2, v = 80000*(highlight.chrom.nums[x] - 1) + cumsum(c(0, chrom.lengths))[highlight.chrom.nums[x]] + highlight.positions[x])
    }
  } else {
    plot(positions[which(chroms == chrom)], lod.vector[which(chroms == chrom)], col = col1, lwd = lwd, xlim = coords, main = title, xlab = paste(chrom, "position"), ylab = "LOD")
    if (!is.na(highlight.markers[1])) {
      highlight.chroms <- sapply(strsplit(highlight.markers, "_"), function (x) x[1])
      highlight.positions <- as.numeric(sapply(strsplit(highlight.markers, "_"), function (x) x[2]))
      if (chrom %in% highlight.chroms) for (x in which(highlight.chroms == chrom)) abline(col = "red", lty = 2, v = highlight.positions[x])
    }
  }
}

trait.lod.plots.single.plot <- function (lod.vector.list, qtls = NA, col1 = "cornflowerblue", col2 = "blue", lwd = 2, ymax = 0, title = NA, lod.axis.rep = length(lod.vector.list)) {
  if (ymax == 0) ymax <- ceiling(max(unlist(lod.vector.list))/10)*10
  plot(1, type = "n", xlim = c(0, sum(chrom.lengths[-17]) + 80000*15), ylim = c(0, (ymax + 20) * length(lod.vector.list)), axes = F, xlab = "chrom", ylab = "LOD", main = title)
  for (x in 1:lod.axis.rep) {
    axis(2, col = NA, col.ticks = 1, at = (ymax + 20) * (length(lod.vector.list) - x) + seq(0, ymax - 20, length.out = 3), labels = seq(0, ymax - 20, length.out = 3), las = 2)
  }
  for (x in 1:16) {
    chrom.lims <- c(seq(0,80000*16, by= 80000) + cumsum(c(0,chrom.lengths))[1:17])[x:(x + 1)] + c(0, -80000)
    axis(1, at = seq(chrom.lims[1], chrom.lims[2], length.out = 3), 
         labels = as.roman(c(NA, x, NA)), 
         col = "#777777", lwd.ticks = 0, lwd = 2, lend = 1)
  }
  for (i in 1:length(lod.vector.list)) {
    lod.plot(lod.vector.list[[i]], names(lod.vector.list[[i]]), y.offset = (ymax + 20) * (length(lod.vector.list) - i), add = T)
  }
  if (is.data.frame(qtls)) {
    for (x in 1:nrow(qtls)) {
      text(qtls$pos[x] + sum(c(0, chrom.lengths)[1:as.numeric(as.roman(substr(qtls$pmarker[x], 4, nchar(qtls$pmarker[x]))))] + 80000) - 80000,
      )
    }
  }
}

max.chrVII.maltose.LODS <- numeric(0)
maltose.lod.vectors <- list()

for (cross in crosses) {
  genotype.file <- paste0("~/data/pacbio RR genomes/input/QTL_data/", list.files("~/data/pacbio RR genomes/input/QTL_data/")[grep(paste0(cross, ".tsv$", collapse = ""), list.files("~/data/pacbio RR genomes/input/QTL_data/"))], collapse = "")
  genos <- fread(genotype.file)
  trait <- grep("maltose", colnames(phenos), ignore.case = T)
  phenos.trait.cross <- phenos[which(phenos$id %in% genos$V1), trait]
  maltose.lod.vectors[[length(maltose.lod.vectors) + 1]] <- lods(phenos.trait.cross, genos)
  names(maltose.lod.vectors[[length(maltose.lod.vectors)]]) <- colnames(genos)[-1]
  max.chrVII.maltose.LODS <- c(max.chrVII.maltose.LODS, max(maltose.lod.vectors[[length(maltose.lod.vectors)]][grep("chrVII_", colnames(genos))]))
}
max.chrVII.maltose.LODS
save(maltose.lod.vectors, file = "~/data/pacbio RR genomes/output/maltose/maltose.lod.vectors.Rdata")
load("~/data/pacbio RR genomes/output/maltose/maltose.lod.vectors.Rdata")

pdf(file = "~/data/pacbio RR genomes/output/maltose/lod.plots.2.pdf", height = 5, width = 3.5)
trait.lod.plots.single.plot(maltose.lod.vectors, ymax = 200)
dev.off()

#Removing the ade- segregants from the last two crosses
ade.marker <- c(which.max(maltose.lod.vectors[[15]]), which.max(maltose.lod.vectors[[16]]))
maltose.lod.vectors.filt <- list()

for (x in 1:2) {
  genotype.file <- paste0("~/data/pacbio RR genomes/input/QTL_data/", list.files("~/data/pacbio RR genomes/input/QTL_data/")[grep(paste0(crosses[14 + x], ".tsv$", collapse = ""), list.files("~/data/pacbio RR genomes/input/QTL_data/"))], collapse = "")
  genos <- fread(genotype.file)
  genos <- as.data.frame(genos)
  genos <- genos[which(genos[, ade.marker[x]] == as.numeric(names(sort(table(genos[,ade.marker[x]]), decreasing = T))[1])),]
  trait <- grep("maltose", colnames(phenos), ignore.case = T)
  phenos.trait.cross <- phenos[which(phenos$id %in% genos$V1), trait]
  maltose.lod.vectors.filt[[length(maltose.lod.vectors.filt) + 1]] <- lods(phenos.trait.cross, genos)
  names(maltose.lod.vectors.filt[[length(maltose.lod.vectors.filt)]]) <- colnames(genos)[-1]
}

save(maltose.lod.vectors.filt, file = "~/data/pacbio RR genomes/output/maltose/maltose.lod.vectors.filt.Rdata")
load("~/data/pacbio RR genomes/output/maltose/maltose.lod.vectors.filt.Rdata")

pdf(file = "~/data/pacbio RR genomes/output/maltose/lod.plots.clib219.pdf", height = 5, width = 3.5)
trait.lod.plots.single.plot(maltose.lod.vectors.filt, ymax = 200) 
dev.off()

#Paraquat

genotype.dir = "~/data/pacbio RR genomes/input/QTL_data/"
paraquat.column <- grep("paraquat", colnames(phenos), ignore.case = T)
paraquat.lod.vector <- list()
for (x in 1:length(crosses)) {
  genotype.file <- paste0(genotype.dir, list.files(genotype.dir)[grep(paste0(crosses[x], ".tsv$", collapse = ""), list.files(genotype.dir))], collapse = "")
  genos <- fread(genotype.file)
  phenos.trait.cross <- phenos[which(phenos$id %in% genos$V1), paraquat.column]
  paraquat.lod.vector[[x]] <- lods(phenos.trait.cross, genos)
  names(paraquat.lod.vector[[x]]) <- colnames(genos)[-1]
}

save(paraquat.lod.vector, file = "~/data/pacbio RR genomes/output/paraquat/paraquat.lod.vector.Rdata")
load("~/data/pacbio RR genomes/output/paraquat/paraquat.lod.vector.Rdata")

pdf(file = "~/data/pacbio RR genomes/output/paraquat/lod.plots.2.pdf", height = 5, width = 3.5)
trait.lod.plots.single.plot(paraquat.lod.vector, ymax = 200)
dev.off()


#Removing the ade- segregants from the last two crosses
ade.marker <- c(which.max(maltose.lod.vectors[[15]]), which.max(maltose.lod.vectors[[16]]))
paraquat.lod.vectors.filt <- list()

for (x in 1:2) {
  genotype.file <- paste0(genotype.dir, list.files(genotype.dir)[grep(paste0(crosses[14 + x], ".tsv$", collapse = ""), list.files(genotype.dir))], collapse = "")
  genos <- fread(genotype.file)
  genos <- as.data.frame(genos)
  genos <- genos[which(genos[, ade.marker[x]] == as.numeric(names(sort(table(genos[,ade.marker[x]]), decreasing = T))[1])),]
  trait <- grep("paraquat", colnames(phenos), ignore.case = T)
  phenos.trait.cross <- phenos[which(phenos$id %in% genos$V1), trait]
  paraquat.lod.vectors.filt[[length(paraquat.lod.vectors.filt) + 1]] <- lods(phenos.trait.cross, genos)
  names(paraquat.lod.vectors.filt[[length(paraquat.lod.vectors.filt)]]) <- colnames(genos)[-1]
}

save(paraquat.lod.vectors.filt, file = "~/data/pacbio RR genomes/output/paraquat/paraquat.lod.vectors.filt.Rdata")
load("~/data/pacbio RR genomes/output/paraquat/paraquat.lod.vectors.filt.Rdata")

pdf(file = "~/data/pacbio RR genomes/output/paraquat/lod.plots.clib219.pdf", height = 5, width = 3.5)
trait.lod.plots.single.plot(paraquat.lod.vectors.filt, ymax = 200) 
dev.off()

#########
#MAL and SGE1 GENE TRACK PLOTS
#########

#MALTOSE

#First, I printed fasta files of chrVII from BIO2 to the chromosome end for every genome.
chrVIIR.seqs <- list()
for (x in 1:16) chrVIIR.seqs[[x]] <- DNAString(dna.sequence(fastas[[x]], gff.files.fams[[x]]$V1[which(gff.files.fams[[x]]$ID == "YGR286C")], gff.files.fams[[x]]$V4[which(gff.files.fams[[x]]$ID == "YGR286C")], length(fastas[[x]][[gff.files.fams[[x]]$V1[which(gff.files.fams[[x]]$ID == "YGR286C")]]])))
names(chrVIIR.seqs) <- paste0("MSY", 24:39, " chrVII-R end")
for (x in 1:16) writeXStringSet(DNAStringSet(chrVIIR.seqs[[x]]), filepath = paste0("~/data/pacbio RR genomes/output/maltose/chrVII_fastas/chrVII-R - MSY", x + 23, ".fasta"), format = "fasta")

#I then imported the sequences into benchling and used blast to manually annotate MAL genes, to produce the chrVII_mal_genes.csv file below. 
#Genes after the last MAL gene are not included here - in some cases there are quite a few, so I don't think they 
#should be included, but it may be nice to give some idea of how close to the telomere end the last MAL gene is.
mal.genes <- read.csv("~/data/pacbio RR genomes/input/maltose/chrVII_mal_genes.csv", header = T, stringsAsFactors = F)

mal.genes.by.strain <- split(mal.genes, mal.genes$Strain..MSY.)

plot.gene <- function (gene.start, gene.end, y = 0, col = "white", label = NA, height = 300, density = NULL, font.size = 1) {
  polygon(c(gene.end, gene.end - (height/2)*(sign(gene.end - gene.start)), gene.start, gene.start, gene.end - (height/2)*(sign(gene.end - gene.start))) ,
          c(y, y+height/2, y + height/2, y - height/2, y - height/2),
          col = col, density = density)
  text(gene.start + 125 * sign(gene.end - gene.start), y, label, adj = c(c(0,1)[sign(gene.end - gene.start)], 0.5), cex = font.size, font = 3)
}

pdf(file = "~/data/pacbio RR genomes/output/maltose/mal.genes.pdf", height = 10, width = 15)
plot(1, type = "n", xlim = c(0, 40000), ylim = c(0, 25000), asp = 1)
sapply(1:length(mal.genes.by.strain), function (y) {
  sapply(1:(nrow(mal.genes.by.strain[[y]]) - 1), function (x) {
    plot.gene(mal.genes.by.strain[[y]]$start[x] - min(mal.genes.by.strain[[y]][,c("start", "end")]),
              mal.genes.by.strain[[y]]$end[x] - min(mal.genes.by.strain[[y]][,c("start", "end")]),
              25000 - 1500*y,
              c("#42a5f5", "yellow", "#800080", "grey")[which(unique(mal.genes$gene) == mal.genes.by.strain[[y]]$gene[x])],
              label = mal.genes.by.strain[[y]]$gene[x],
              height = 700,
              font.size = 0.7,
              density = if (identical(mal.genes.by.strain[[y]]$pseudogene_.non.exhaustive_list.[x], T)) 40 else NULL)
  })
})
dev.off()

#Extra MAL genes on chromosome XI
chrXIR.seqs <- list()
for (x in 1:16) chrXIR.seqs[[x]] <- DNAString(dna.sequence(fastas[[x]], gff.files.fams[[x]]$V1[which(gff.files.fams[[x]]$ID == "YKR101W")], gff.files.fams[[x]]$V4[which(gff.files.fams[[x]]$ID == "YKR101W")], length(fastas[[x]][[gff.files.fams[[x]]$V1[which(gff.files.fams[[x]]$ID == "YKR101W")]]])))
names(chrXIR.seqs) <- paste0("MSY", 24:39, " chrXI-R end")
for (x in 1:16) writeXStringSet(DNAStringSet(chrXIR.seqs[[x]]), filepath = paste0("~/data/pacbio RR genomes/output/maltose/chrXI_fastas/chrXI-R - MSY", x + 23, ".fasta"), format = "fasta")
#Annotated using benchling and blast, and then made graphs in illustrator from the manual annotation directly, without using R to make gene plots as above.

#SGE1. Annotated using Benchling and BLAST.
other.paraquat.genes <- read.csv("~/data/pacbio RR genomes/input/paraquat/chrVII_and_XI_paraquat_genes.csv", stringsAsFactors = F, header = T)
other.paraquat.genes.by.strain <- split(other.paraquat.genes, other.paraquat.genes$Strain..MSY.)

pdf(file = "~/data/pacbio RR genomes/output/paraquat/other.sge1.genes.pdf", height = 10, width = 15)
plot(1, type = "n", xlim = c(0, 40000), ylim = c(0, 25000), asp = 1)
sapply(1:length(other.paraquat.genes.by.strain), function (y) {
  sapply(1:(nrow(other.paraquat.genes.by.strain[[y]]) - 1), function (x) {
    plot.gene(other.paraquat.genes.by.strain[[y]]$start[x] - min(other.paraquat.genes.by.strain[[y]][,c("start", "end")]),
              other.paraquat.genes.by.strain[[y]]$end[x] - min(other.paraquat.genes.by.strain[[y]][,c("start", "end")]),
              25000 - 1500*y,
              c("pink", "lavender", "white", "#42a5f5", "orange", "yellow", "#800080")[which(unique(other.paraquat.genes$gene) == other.paraquat.genes.by.strain[[y]]$gene[x])],
              label = other.paraquat.genes.by.strain[[y]]$gene[x],
              height = 700,
              font.size = 0.7,
              density = if (identical(other.paraquat.genes.by.strain[[y]]$pseudogene_.non.exhaustive_list.[x], T)) 40 else NULL)
  })
})
dev.off()


########
#SCATTERPLOTS OF SEGREGANT PHENOTYPES
########

segregant.t.test <- function (cross, phenotype, test.chrom, test.site, dependent.chroms = character(0), dependent.sites = numeric(0), filt.chrom = NA, filt.site = NA, filt.geno = NA, figure = F) {
  if (grepl("MSY", cross)) {
    cross <- c("375", "A", "376", "B", "377", "393", "381", "3008", "2999", "3000", "3001", "3049", "3003", "3004", "3043", "3028")[as.numeric(substr(cross, 4, 5)) - 23]
  }
  genos <- fread(paste0("~/data/pacbio RR genomes/input/QTL_data/genotype_", cross, ".tsv", collapse = ""))
  if (!is.na(filt.chrom)) {
    filt.marker <- names(which.min(abs(sapply(grep(paste0(filt.chrom, "_"), colnames(genos), value = T), function (x) as.numeric(strsplit(x, "_")[[1]][2])) - filt.site)))
    genos <- genos[which(as.data.frame(genos)[, filt.marker] == filt.geno),]
  }
  geno.test.column <- names(which.min(abs(sapply(grep(paste0(test.chrom, "_"), colnames(genos), value = T), function (x) as.numeric(strsplit(x, "_")[[1]][2])) - test.site)))
  if (length(dependent.chroms) > 0) {
    geno.dep.columns <- sapply(1:length(dependent.chroms), function (y) names(which.min(abs(sapply(grep(paste0(dependent.chroms[y], "_"), colnames(genos), value = T), function (x) as.numeric(strsplit(x, "_")[[1]][2])) - dependent.sites[y]))))
  } else geno.dep.columns <- numeric(0)
  trait <- grep(phenotype, colnames(phenos), ignore.case = T)
  phenos.trait.cross <- phenos[which(phenos$id %in% genos$V1), trait]
  phenos.by.genos <- split(phenos.trait.cross, as.data.frame(genos)[,c(geno.test.column, geno.dep.columns)])
  t.test.results <- sapply(1:(length(phenos.by.genos)/2), function (x) t.test(phenos.by.genos[[x * 2 - 1]], phenos.by.genos[[x * 2]]))
  if (figure == T) {
    stripchart(phenos.by.genos, vertical = T, method = "jitter", jitter = .2, pch = 1, col = "#00000055", ylab = paste0("Growth on ", phenotype), xlab = paste(c("Genotype at ", paste(c(geno.test.column, geno.dep.columns), collapse = ", ")), collapse = ""), axes = F, main = cross)
    axis(2)
    axis(1, at = 1:length(phenos.by.genos), names(phenos.by.genos), las = 2, cex.axis = 0.7)
    for (x in 1:length(phenos.by.genos)) lines(c(x - 0.3, x + 0.3), c(mean(phenos.by.genos[[x]]), mean(phenos.by.genos[[x]])), lwd = 2, col = "blue")
  }
  t.test.results
}

pdf("~/data/pacbio RR genomes/output/paraquat/segregant scatterplots.pdf", height = 5.5, width = 8.5)
segregant.t.test("MSY25_MSY26", "paraquat", "chrXI", 666499, figure = T) #p-value 0.0018
segregant.t.test("MSY26_MSY27", "paraquat", "chrXI", 666499, "chrXVI", 946477, figure = T) #Neither is significant: p-value 0.8 and 0.98
segregant.t.test("MSY35_MSY36", "paraquat", "chrXI", 666499, "chrXVI", 946477, figure = T) #Neither is significant: p-value 0.55 and 0.27
segregant.t.test("MSY36_MSY37", "paraquat", "chrXI", 666499, "chrXVI", 946477, figure = T) #Both are significant: p-value 3.8e-5 and 6.3e-27
segregant.t.test("MSY38_MSY39", "paraquat", "chrVII", 1083575, "chrXVI", 946477) #Both are significant: p-value 0.00089 and 4.6e-7. These are not filtered to remove the ade- colonies.
segregant.t.test("MSY39_MSY24", "paraquat", "chrVII", 1083575, "chrXVI", 946477) #Both are significant: p-value 1.1e-10 and 3.7e-50. These are not filtered to remove the ade- colonies.

segregant.t.test("MSY38_MSY39", "paraquat", "chrVII", 1083575, "chrXVI", 946477, filt.chrom = "chrXV", filt.site = 565236, filt.geno = 1, figure = T) #Still significant after filtering. p-value 0.0030 and 1.1e-7
segregant.t.test("MSY39_MSY24", "paraquat", "chrVII", 1083575, "chrXVI", 946477, filt.chrom = "chrXV", filt.site = 565236, filt.geno = 2, figure = T) #Still significant after filtering. p-value 4.8e-10 and 2.1e-51
dev.off()

pdf("~/data/pacbio RR genomes/output/maltose/segregant scatterplots.pdf", height = 5.5, width = 8.5)
segregant.t.test("MSY25_MSY26", "maltose", "chrXI", 666499, "chrVII", 1083575, figure = T)
segregant.t.test("MSY26_MSY27", "maltose", "chrXI", 666499, "chrVII", 1083575, figure = T)
segregant.t.test("MSY28_MSY29", "maltose", "chrXI", 666499, "chrVII", 1083575, figure = T)
segregant.t.test("MSY29_MSY30", "maltose", "chrXI", 666499, "chrVII", 1083575, figure = T)
dev.off()

pdf("~/data/pacbio RR genomes/output/sucrose/segregant scatterplots - 2.pdf", height = 5.5, width = 8.5)
segregant.t.test("MSY36_MSY37", "sucrose", "chrIV", 9228, "chrIX", 33216, figure = T)
segregant.t.test("MSY35_MSY36", "sucrose", "chrIV", 9228, figure = T)
segregant.t.test("MSY27_MSY28", "sucrose", "chrX", 734496, figure = T)
segregant.t.test("MSY36_MSY37", "raffinose", "chrIV", 9228, "chrIX", 33216, figure = T)
segregant.t.test("MSY35_MSY36", "raffinose", "chrIV", 9228, figure = T)
segregant.t.test("MSY27_MSY28", "raffinose", "chrX", 734496, figure = T)
dev.off()

###########
#MALR Percent ID PLOTS
###########

dnapid <- read.csv("~/data/pacbio RR genomes/input/maltose/chrVII-R MALR dna pid 221024.csv", row.names = 1, header = T, stringsAsFactors = F) #Generated by CLUSTAL omega from MALR sequences pulled from the benchling annotation mentioned above.

#DNA percent ID
pdf(file = "~/data/pacbio RR genomes/output/maltose/chrVII-R MALR dna pid 2201024.pdf", height = 11, width = 8.5)
pheatmap(as.matrix(dnapid), border_color = NA, cellwidth = 10, cellheight = 10, breaks = c(seq(min(dnapid), 99, length.out = 90), seq(99, 100, length.out = 12)[-1]))
dev.off()

#Within a given strain's chrVII MAL cluster, how similar are the MALRs?
strains.w.chrvii.malrs.dna <- unique(sapply(strsplit(colnames(dnapid), "\\."), function (x) strsplit(x[1], "MALR_")[[1]][2]))
strains.w.chrvii.malrs.dna <- strains.w.chrvii.malrs.dna[-which(is.na(strains.w.chrvii.malrs.dna))]
sapply(strains.w.chrvii.malrs.dna, function (y) dnapid[grep(y, colnames(dnapid)), grep(y, colnames(dnapid))])

dnapid.nodiag <- dnapid
diag(dnapid.nodiag) <- NA
median(unlist(sapply(strains.w.chrvii.malrs.dna, function (y) unlist(dnapid.nodiag[grep(y, colnames(dnapid)), grep(y, colnames(dnapid))]))), na.rm = T)

###########
#Chr XI MALR chimerism graph
###########

chrximalrs <- readAAMultipleAlignment("~/data/pacbio RR genomes/input/maltose/chrxi_malrs_alignment.clustal_num", format = "clustal") #Alignment from clustal omega 2.1 generated from benchling annotated genes
mal23 <- strsplit(as.character(chrximalrs@unmasked$MAL23), "")[[1]]
malr.msy26 <- strsplit(as.character(chrximalrs@unmasked$RM_chrXI), "")[[1]]
malr.msy29 <- strsplit(as.character(chrximalrs@unmasked$MSY29_chrXI), "")[[1]]
mal43 <- strsplit(as.character(chrximalrs@unmasked$MAL43), "")[[1]]

plot(type = "n", xlim = c(0, length(mal23)), ylim = c(-1, 5), axes = F, 1, xlab = "", ylab = "")

for (x in 1:(length(mal23) - 1)) {
  if ((mal23[x] == mal43[x]) && (malr.msy26[x] == mal23[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "white", border = NA)
  } else if (malr.msy26[x] == mal23[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "green", border = NA)
  } else if (malr.msy26[x] == mal43[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "purple", border = NA)
  } else {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "orange", border = NA)
  }
}

for (x in 1:(length(mal23) - 1)) {
  if ((mal23[x] == mal43[x]) && (malr.msy29[x] == mal23[x])) {
    polygon(c(x, x + 1, x + 1, x), c(1.2, 1.2, 1.8, 1.8), col = "white", border = NA)
  } else if (malr.msy29[x] == mal23[x]) {
    polygon(c(x, x + 1, x + 1, x), c(1.2, 1.2, 1.8, 1.8), col = "green", border = NA)
  } else if (malr.msy29[x] == mal43[x]) {
    polygon(c(x, x + 1, x + 1, x), c(1.2, 1.2, 1.8, 1.8), col = "purple", border = NA)
  } else {
    polygon(c(x, x + 1, x + 1, x), c(1.2, 1.2, 1.8, 1.8), col = "orange", border = NA)
  }
}

for (x in 1:(length(mal23) - 1)) {
  if (mal23[x] == mal43[x]) {
    polygon(c(x, x + 1, x + 1, x), c(-.2, -.2, -.8, -.8), col = "white", border = NA)
  } else {
    polygon(c(x, x + 1, x + 1, x), c(-.2, -.2, -.8, -.8), col = "green", border = NA)
  }
}

for (x in 1:(length(mal23) - 1)) {
  if (mal23[x] == mal43[x]) {
    polygon(c(x, x + 1, x + 1, x), c(2.2, 2.2, 2.8, 2.8), col = "white", border = NA)
  } else {
    polygon(c(x, x + 1, x + 1, x), c(2.2, 2.2, 2.8, 2.8), col = "purple", border = NA)
  }
}

sapply(0:3, function (x) polygon(c(1, length(mal23), length(mal23), 1), c(-.2 + x, -.2 + x, -.8 + x, -.8 + x), col = NA))

axis(2, at = seq(-.5, 2.5), c("MAL23", "MALR_RM.chrXI", "MALR_CLIB413.chrXI", "MAL43"), las = 2, tick = F)

###########
#SUC gene locations
###########

cen.positions <- read.csv("~/data/pacbio RR genomes/input/centromere_positions.csv", stringsAsFactors = F) #Manually curated from SGD.

suc.gene.locations <- do.call("rbind", lapply(1:16, function (x) gff.files.fams[[x]][grep("YIL162W", gff.files.fams[[x]]$V9),-9]))
suc.gene.copy.number <- sapply(1:16, function (x) length(grep("YIL162W", gff.files.fams[[x]]$V9)))
suc.gene.locations <- cbind.data.frame(suc.gene.locations, strain = rep(strains, suc.gene.copy.number))
suc.gene.locations <- suc.gene.locations[order(suc.gene.locations$CHROM),]

plot(1, type = "n", 
     xlim = c(0, max(chrom.lengths)), 
     ylim = c(0, 34), 
     axes = F, ylab = "", # ylab = "Chromosome",
     xlab = "Chromosome position (kb)", )

for (x in 1:16) {
  polygon(c(0, cen.positions$left[x], cen.positions$left[x], 0),
          c(33-x*2, 33 - x*2, 32 - x*2, 32 - x*2),
          col = "#FFFFEE")
  polygon(c(cen.positions$right[x], chrom.lengths[x], chrom.lengths[x], cen.positions$right[x]),
          c(33-x*2, 33 - x*2, 32 - x*2, 32 - x*2),
          col = "#FFFFEE")
}

axis(side = 1, at = seq(0, 1500000, by = 500000), labels = c(0, 500, 1000, 1500), col = NA, col.ticks = "#000000")


for (x in 1:nrow(suc.gene.locations)) lines(
  col = "blue",
  c(suc.gene.locations$V4[x], suc.gene.locations$V4[x]), 
  c(33 - 2 * as.numeric(as.roman(substr(suc.gene.locations$CHROM[x], 4, nchar(as.character(suc.gene.locations$CHROM[x]))))),
    32 - 2 * as.numeric(as.roman(substr(suc.gene.locations$CHROM[x], 4, nchar(as.character(suc.gene.locations$CHROM[x])))))))
#I marked chromosome ends with the strains having telomeric SUCs in illustrator, and consolidated the SUC2 chrIX loci into a single mark.
#I also changed to rounded rectangles for the chromosome arms.

for (x in 1:16) {
  #  lines(col = "#AAAAAA", rep(.5 * (cen.positions[,2] + cen.positions[,3])[x], 2), c(33 - 2 * x, 32 - 2 * x))
  text(.5 * (cen.positions[,2] + cen.positions[,3])[x] + 10000 * (if (x %in% c(1, 6, 9, 14)) -1 else 1), 32.5 - 2 * x, paste0("chr", as.roman(x)), adj = c((if (x %in% c(1, 6, 9, 14)) 1 else 0), 0.5))
}

#######
#Genetic differences between SUC genes
#######

sucgenes <- DNAStringSet(unlist(sapply(names(fastas), function (x) {
  sapply(grep("YIL162W", gff.files[[grep(x, names(gff.files), ignore.case = T)]]$V9), function (y) {
    DNAString(dna.sequence.from.gff(x, y))
  })
})))
names(sucgenes) <- unlist(sapply(1:16, function (x) paste0(names(fastas)[x], "_", gff.files[[x]]$CHROM[grep("YIL162W", gff.files[[x]]$V9)], "_", gff.files[[x]]$V4[grep("YIL162W", gff.files[[x]]$V9)])))
writeXStringSet(sucgenes, filepath = "~/data/pacbio RR genomes/output/sucrose/suc.fasta", format = "fasta")
sucgenes <- readDNAStringSet("~/data/pacbio RR genomes/output/sucrose/suc_wpar.fasta") #I added paradoxus SUC2 (obtained from SGD) to the fasta file.

for (x in 2:length(sucgenes)) {
  if (grepl("chrIX", names(sucgenes)[x])) {
    names(sucgenes)[x] <- paste0("SUC2_", strains[as.numeric(substr(names(sucgenes)[x], 4, 5)) - 23])
  } else {
    names(sucgenes)[x] <- paste0("telSUC_", strains[as.numeric(substr(names(sucgenes)[x], 4, 5)) - 23], ".", strsplit(names(sucgenes)[x], "_")[[1]][2], "-", c("L", "R")[sign(100000 - as.numeric(strsplit(names(sucgenes)[x], "_")[[1]][3]))])
  }
}
names(sucgenes)[13] <- "telSUC_YJM454.chrIX-L" #There's one SUC gene on chromosome 9 that's actually a telomeric SUC. All other chr IX SUCs are at the SUC2 locus.

#Percent ID heatmap
suc.dna.pid <- read.csv("~/data/pacbio RR genomes/input/sucrose/suc dna pid.csv", row.names = 1, header = T, stringsAsFactors = F) #Generated by CLUSTAL omega2.1

pdf(file = "~/data/pacbio RR genomes/output/sucrose/suc.dna.pid.pdf", height = 11, width = 8.5)
pheatmap(as.matrix(suc.dna.pid), border_color = NA, cellwidth = 10, cellheight = 10, breaks = c(seq(min(dnapid), 99, length.out = 90), seq(99, 100, length.out = 12)[-1]))
dev.off()

#Takeaways from the clustal percent ID table:
#Overall identity between the dominant class of telomeric SUC genes and the reference SUC2 allele is 93%.
median(sort(suc.dna.pid[,'SUC2_BY'])[1:13])

#Divergence among the dominant class of subtelomeric SUC genes is low; percent ID is all above 99.5%
min(suc.dna.pid[grep("telSUC", rownames(suc.dna.pid))[1:11], grep("telSUC", rownames(suc.dna.pid))[1:11]])

suc2.by <- strsplit(as.character(sucgenes$SUC2_BY), "")[[1]]
suc2.paradoxus <- strsplit(as.character(sucgenes$paradoxus), "")[[1]]
suc2.273614 <- strsplit(as.character(sucgenes$SUC2_273614), "")[[1]]
telsuc.yps1009 <- strsplit(as.character(sucgenes$`telSUC_YPS1009.chrVII-R`), "")[[1]]
telsuc.273614 <- strsplit(as.character(sucgenes$`telSUC_273614.chrIV-L`), "")[[1]]

#Telomeric SUC (SUC5 from 273614) vs SUC2
plot(type = "n", xlim = c(0, length(sucgenes$SUC2_BY)), ylim = c(-1, 5), axes = F, 1, xlab = "", ylab = "")
for (x in 1:length(sucgenes$SUC2_BY)) {
  if (suc2.by[x] == telsuc.273614[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "white", border = NA)
  } else if (suc2.by[x] != telsuc.273614[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "black", border = NA)
  }
}
polygon(c(1, 1599, 1599, 1), c(.2, .2, .8, .8))

#Levels of sequence identity in regions of concentrated SNPs
1 - length(which(suc2.by[1:432] != telsuc.273614[1:432]))/432 #88.0% identity from 1-432
1 - length(which(suc2.by[535:718] != telsuc.273614[535:718]))/(718 - 535 + 1) #91.8% identity from 535 to 718
1 - length(which(suc2.by[1053:1292] != telsuc.273614[1053:1292]))/(1292 - 1053 + 1) #81.7% identity from 1053 to 1292

#Paradoxus vs. cerevisiae SUC2
plot(type = "n", xlim = c(0, length(sucgenes$SUC2_BY)), ylim = c(-1, 5), axes = F, 1, xlab = "", ylab = "")
for (x in 1:length(sucgenes$SUC2_BY)) {
  if (suc2.by[x] == suc2.paradoxus[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "white", border = NA)
  } else if (suc2.by[x] != suc2.paradoxus[x]) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "black", border = NA)
  }
}
polygon(c(1, 1599, 1599, 1), c(.2, .2, .8, .8))

#Levels of sequence identity in the same regions as above
1 - length(which(suc2.by[1:432] != suc2.paradoxus[1:432]))/432 #92.6% identity from 1-432
1 - length(which(suc2.by[535:718] != suc2.paradoxus[535:718]))/(718 - 535 + 1) #89.7% identity from 535 to 718
1 - length(which(suc2.by[1053:1292] != suc2.paradoxus[1053:1292]))/(1292 - 1053 + 1) #86.3% identity from 1053 to 1292

#Comparing all three at once
plot(type = "n", xlim = c(0, length(sucgenes$SUC2_BY)), ylim = c(-1, 5), axes = F, 1, xlab = "", ylab = "")
for (x in 1:length(sucgenes$SUC2_BY)) {
  if ((suc2.by[x] == suc2.paradoxus[x]) && (suc2.by[x] == telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "white", border = NA)
  } else if ((suc2.by[x] != suc2.paradoxus[x]) && (suc2.by[x] == telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "green", border = NA)
  } else if ((suc2.by[x] == suc2.paradoxus[x]) && (suc2.by[x] != telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "purple", border = NA)
  } else if ((suc2.by[x] != suc2.paradoxus[x]) && (suc2.paradoxus[x] == telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "orange", border = NA)
  } else if (((suc2.by[x] != suc2.paradoxus[x]) && (suc2.paradoxus[x] != telsuc.273614[x])) && (suc2.by[x] != telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "black", border = NA)
  }
}
polygon(c(1, 1599, 1599, 1), c(.2, .2, .8, .8))

#Comparing the second type of subtelomeric SUC gene, found in YPS1009 and at the SUC2 locus for 273614, to BY_SUC2 and the telomeric SUC in 273614.
plot(type = "n", xlim = c(0, length(sucgenes$SUC2_BY)), ylim = c(-1, 5), axes = F, 1, xlab = "", ylab = "")
for (x in 1:length(sucgenes$SUC2_BY)) {
  if ((suc2.by[x] == telsuc.yps1009[x]) && (suc2.by[x] == telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "white", border = NA)
  } else if ((suc2.by[x] != telsuc.yps1009[x]) && (telsuc.yps1009[x] != telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "orange", border = NA)
  } else if ((suc2.by[x] == telsuc.yps1009[x]) && (suc2.by[x] != telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "purple", border = NA)
  } else if ((suc2.by[x] != telsuc.yps1009[x]) && (telsuc.yps1009[x] == telsuc.273614[x])) {
    polygon(c(x, x + 1, x + 1, x), c(0.2, 0.2, .8, .8), col = "green", border = NA)
  } 
}
polygon(c(1, 1599, 1599, 1), c(.2, .2, .8, .8))

#############
#MALTOSE GROWTH CURVES
#############

plot.growth.curve <- function (strains.to.plot, growth.data, neg.ctrl.data, key.vec, colors = c("black", "deepskyblue", "grey70", "gold"), smooth = 0.000005, time.interval = 1, axis.interval = if (time.interval == 1) 10 else 24, max.timepoint = NA, title = "Growth curves", lwd = 2, yrange = NA, add = F) {
  if (is.na(max.timepoint)) max.timepoint <- ncol(growth.data) - 1
  data.lines <- list()
  for (i in (1:length(strains.to.plot))) data.lines[[i]] <- grep(strains.to.plot[i], key.vec)
  if (length(yrange) == 1) yrange <- range(unlist(growth.data[unlist(data.lines), -1]))
  if (add == F) plot(axes = F, ylab = "OD600", xlab = if (time.interval == 1) "Timepoints" else "Time (hrs)", ylim = yrange, xlim = c(0, max.timepoint), type = "n", 1, main = title)
  for (y in 1:length(data.lines)) {
    for (x in 1:length(data.lines[[y]])) {
      if (length(neg.ctrl.data) > 1) {
        lines(smooth.spline(as.numeric(growth.data[data.lines[[y]][x], -1]) - as.numeric(neg.ctrl.data[data.lines[[y]][x], -1]), lambda = smooth), col = colors[y], lwd = lwd)
      } else if (neg.ctrl.data == "t0")
        lines(smooth.spline(as.numeric(growth.data[data.lines[[y]][x], -1]) - as.numeric(growth.data[data.lines[[y]][x], 2]), lambda = smooth), col = colors[y], lwd = lwd)
      else lines(smooth.spline(as.numeric(growth.data[data.lines[[y]][x], -1]), lambda = smooth), col = colors[y], lwd = lwd)
    }
  }
  axis(2)
  axis(1, seq(0, 6000, by = axis.interval*60)/time.interval, labels = seq(0, 100, by = axis.interval))
}

glucose0418 <- read.csv("~/data/pacbio RR genomes/input/220418 - MAL63 plate reader/plate reader data/glucose.CSV", stringsAsFactors = F)
nosugar0418 <- read.csv("~/data/pacbio RR genomes/input/220418 - MAL63 plate reader/plate reader data/no sugar.CSV", stringsAsFactors = F)
maltose0418 <- read.csv("~/data/pacbio RR genomes/input/220418 - MAL63 plate reader/plate reader data/maltose.CSV", stringsAsFactors = F)

key0418 <- read.csv("~/data/pacbio RR genomes/input/220418 - MAL63 plate reader/plate arrangement.csv", stringsAsFactors = F, header = F)
key0418$V11[1] <- "empty" #There were only 4 colonies on the MSY30 + MSp183 transformation, so I left one of its 5 wells blank.
key.vector0418 <- c(t(as.matrix(key0418)))

#Full time-course AUC values
auc.ratios0418full <- (apply(maltose0418[,-1], 1, function (x) sum(as.numeric(x))) - apply(nosugar0418[,-1], 1, function (x) sum(as.numeric(x))))/ (apply(glucose0418[,-1], 1, function (x) sum(as.numeric(x))) - apply(nosugar0418[,-1], 1, function (x) sum(as.numeric(x))))

stripchart(split(auc.ratios0418full, key.vector0418)[1:16], method = "jitter", vertical = T, las = 2, pch = 1, xlab = "Strain with or without MAL63", ylab = "Maltose growth (AUC_maltose/AUC_glucose)")
sapply(1:16, function (x) lines(c(x-.3, x + .3), c(mean(split(auc.ratios0418full, key.vector0418)[[x]]), mean(split(auc.ratios0418full, key.vector0418)[[x]])), col = "blue"))

par(mfrow = c(2, 3))
sapply(c(25, 30, 32, 34, 35, 39), function (x) plot.growth.curve(c(paste0(x, "\\(\\-\\)"), paste0(x, "\\(\\+\\)")), maltose0418, nosugar0418, key.vec = key.vector0418, time.interval = 14.017632, smooth = 0.0001, yrange = c(-0.1, 1.5)))
par(mfrow = c(1,1))

#Glucose curves
par(mfrow = c(2, 3))
sapply(c(25, 30, 32, 34, 35, 39), function (x) plot.growth.curve(c(paste0(x, "\\(\\-\\)"), paste0(x, "\\(\\+\\)")), glucose0418, nosugar0418, key.vec = key.vector0418, time.interval = 14.017632, smooth = 0.0001, yrange = c(-0.1, 1.5)))
par(mfrow = c(1,1))

samples.to.analyze <- intersect(which(nchar(key.vector0418) < 6) , grep("25|25|30|32|34|35|39", key.vector0418))
model.full <- aov(auc.ratios0418full[samples.to.analyze] ~ key.vector0418[samples.to.analyze])
summary(model.full)
TukeyHSD(model.full)

###########
#PARAQUAT GROWTH CURVES
###########

sge1key.plate1 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/base plate arrangement - paraquat 0518.csv", stringsAsFactors = F, header = F)

sge1key.plate1.vector <- c(t(as.matrix(sge1key.plate1)))

noparaquat <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/no paraquat 220518.csv")
paraquat.plate1 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p1 220518.csv")
paraquat.plate2 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p2 220518.csv")
paraquat.plate3 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p3 220518.csv")
paraquat.plate4 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p4 220518.csv")
paraquat.plate5 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p5 220518.csv")
paraquat.plate6 <- read.csv("~/data/pacbio RR genomes/input/220518 - SGE1 across multiple plates/p6 220518.csv")

#Taking the Nth row of each plate (plus the first and last row of water wells) so that it reconstitutes the format of the original plate.
paraquat.row2 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate1[13:24,], paraquat.plate6[13:24,], paraquat.plate5[13:24,], paraquat.plate4[13:24,], paraquat.plate3[13:24,], paraquat.plate2[13:24,], paraquat.plate1[85:96,])
paraquat.row3 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate2[25:36,], paraquat.plate1[25:36,], paraquat.plate6[25:36,], paraquat.plate5[25:36,], paraquat.plate4[25:36,], paraquat.plate3[25:36,], paraquat.plate1[85:96,])
paraquat.row4 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate3[37:48,], paraquat.plate2[37:48,], paraquat.plate1[37:48,], paraquat.plate6[37:48,], paraquat.plate5[37:48,], paraquat.plate4[37:48,], paraquat.plate1[85:96,])
paraquat.row5 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate4[49:60,], paraquat.plate3[49:60,], paraquat.plate2[49:60,], paraquat.plate1[49:60,], paraquat.plate6[49:60,], paraquat.plate5[49:60,], paraquat.plate1[85:96,])
paraquat.row6 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate5[61:72,], paraquat.plate4[61:72,], paraquat.plate3[61:72,], paraquat.plate2[61:72,], paraquat.plate1[61:72,], paraquat.plate6[61:72,], paraquat.plate1[85:96,])
paraquat.row7 <- rbind.data.frame(paraquat.plate1[1:12,], paraquat.plate6[73:84,], paraquat.plate5[73:84,], paraquat.plate4[73:84,], paraquat.plate3[73:84,], paraquat.plate2[73:84,], paraquat.plate1[73:84,], paraquat.plate1[85:96,])

aucs.noparaquat.96hrs <- apply(noparaquat[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row2.96hrs <- apply(paraquat.row2[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row3.96hrs <- apply(paraquat.row3[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row4.96hrs <- apply(paraquat.row4[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row5.96hrs <- apply(paraquat.row5[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row6.96hrs <- apply(paraquat.row6[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))
aucs.row7.96hrs <- apply(paraquat.row7[,2:181], 1, function (x) sum(as.numeric(x)) - as.numeric(x[1])*length(x))

stripchart(split(aucs.row7.96hrs/aucs.noparaquat.96hrs, sge1key.plate1.vector)[-5], vertical = T, method = "jitter", pch = 1, ylim = c(0, 1.2), ylab = "Ratio of AUC_paraquat to AUC_no_paraquat", xlab = "sge1D + SGE1 plasmid")
sapply(c(1:6), function (x) lines(c(x - 0.3, x + 0.3), rep(mean(split(aucs.row7.96hrs/aucs.noparaquat.96hrs, sge1key.plate1.vector)[-5][[x]]), 2), col = "blue"))

model.96hrs <- aov((aucs.row7.96hrs/aucs.noparaquat.96hrs)[-which(sge1key.plate1.vector == "h2o")] ~ sge1key.plate1.vector[-which(sge1key.plate1.vector == "h2o")])
summary(model.96hrs)
TukeyHSD(model.96hrs)

par(mfrow = c(2, 3))
for (x in c("25", "28", "32", "33", "L449M")) plot.growth.curve(c("MSp8", x), paraquat.row7, "t0", key.vec = sge1key.plate1.vector, yrange = c(0, 1.5), max.timepoint = 96*60/32.113744, time.interval = 32.113744, smooth = 0.0001, colors = c("#44444444", "orange"), title = x)
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
for (x in c("25", "28", "32", "33", "L449M")) plot.growth.curve(c("MSp8", x), noparaquat, "t0", key.vec = sge1key.plate1.vector, yrange = c(0, 1.75), max.timepoint = 96*60/32.113744, time.interval = 32.113744, smooth = 0.0001, colors = c("#44444444", "orange"), title = x)
par(mfrow = c(1,1))
