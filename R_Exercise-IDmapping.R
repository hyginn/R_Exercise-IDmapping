# R_Exercise-IDmapping.R
#
# Purpose: Introduction to the intricacies of mapping biological identifiers
#          while creating a mapping file from ENSP IDs to gene symbols.
#
# Version: 1.0
#
# Date:    2017  02  01
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First code
#
# TODO:
#
#
# == HOW TO WORK WITH THIS FILE ================================================
#
#  Go through this script line by line to read and understand the
#  code. Execute code by typing <cmd><enter>. When nothing is
#  selected, that will execute the current line and move the cursor to
#  the next line. You can also select more than one line, e.g. to
#  execute a block of code, or less than one line, e.g. to execute
#  only the core of a nested expression.
#
#  The project directory contains a file "myScript.R". This file is not
#  under version control in the repository and will therefore not create
#  conflicts on updates if you edit it. However, do not _commit_ this file.
#  Use "myScript.R" to write your own code examples, and practice and play.
#  Especially play.
#
#  DO NOT simply source() this whole file!
#
#  If there are portions you don't understand, use R's help system,
#  Google for an answer, or ask me. Don't continue if you don't
#  understand what's going on. That's not how it works ...
#
#
# ==============================================================================
#                  INTRODUCTION
# ==============================================================================


# ID mapping is our daily bread and the bane of our existence. Not only is the
# relationship between genome, gene, transcript, nascent- and native protein
# complex, both the databases and their identifiers  are moving targets because
# they are constantly being updated, and since our data is often quite large,
# there are performance issues to consider - in particular it is usually NOT
# feasible to map IDs over the Internet via the APIs of the database providers.

# The task we will address here is to map from the ENSP IDs that the STRING
# database uses to HGNC gene symbols. Since the human proteins come from the
# human genome, there should be a good correspondence between proteins and
# symbols, right?

# Right?


# Well ...

# ==============================================================================
#        PART ONE: ENSP TO UNIPROT
# ==============================================================================

# The STRING database provides a mapping file from its identifiers to UniProt
# IDs, for many model organisms.
# http://string-db.org/mapping_files/uniprot_mappings/ I have downloaded it and
# put a gzip compressed version in the data directory. Note that R is smart
# enough to be able to read a gzipped file directly into a dataframe - I don't
# have to uncompress the file first.

tsv <- read.delim("data/9606_reviewed_uniprot_2_string.04_2015.tsv.gz",
                  header = FALSE,
                  skip = 1,
                  stringsAsFactors = FALSE)

head(tsv)
IDs <- data.frame(UniProt = gsub("\\|.+$", "", tsv$V2),
                  ENSP = tsv$V3,
                  percID = as.numeric(tsv$V4),
                  stringsAsFactors = FALSE)

head(IDs)
# First, check whether we have entries for each row in all columns
any(is.na(IDs$UniProt) | IDs$UniProt == "")   # FALSE
any(is.na(IDs$ENSP) | IDs$ENSP == "")         # FALSE


# Next, see if there are redundancies
nID <- nrow(IDs)
nID - length(unique(IDs$UniProt))        # 1658
nID - length(unique(IDs$ENSP))           # 1389
# 1,389 of 20,499 records have duplicate ENSP IDs. That's bad, because we need
# to resolve duplicates to be able to map uniquely between identifiers. What do
# these duplicates look like?

dupENS <- which(duplicated(IDs$ENSP))
cat(IDs[IDs$ENSP == IDs$ENSP[dupENS[1]], "UniProt"])
# In this case, the 5 rows with ENSP00000366005 map to UniProt IDs P30443,
# P04439, P13746, P16188, and P30455. If we use the UniProt ID mapping service
# http://www.uniprot.org/uploadlists/ to and convert these UniProtKB AC/ID
# identifiers to "Gene name" (i.e. gene symbols), we see that they are all HLA-A
# - i.e. different variants of the Human histocompatibility antigen.

IDs[IDs$ENSP == IDs$ENSP[dupENS[5]], ]
# The 8 ENSP00000399168 IDs all map to HLA-B

IDs[IDs$ENSP == IDs$ENSP[dupENS[12]], ]
# These 2 ENSP00000353099 IDs map to HLA-DRB1.

# Up to now, it seems that these are different genes of the immune system
# families, and each of these groups has one gene symbol with 100% identity to
# the ENSP ID. However this is not always the case. For example

IDs[IDs$ENSP == IDs$ENSP[dupENS[100]], ]
# ... finds four proteins for ENSP00000401363 which map to DAZ1, DAZ2, DAZ3 and
# DAZ4, a Y-chromosome protein implicated in azoospermia. DAZ1 and 4 appear to
# have the same sequence but DAZ1 is located on Y:NC_000024.10
# (23129355..23199117, complement) whereas DAZ4 is NC_000024.10
# (24833807..24907040). These are recent duplicates.

# Let's check whether all duplicates have at least one UniProt ID with 100%
# sequence ID to the ENSP entry. First, we make a unique list of duplicated ENSP
# entries.

uDupENSP <- unique(IDs$ENSP[dupENS])

N <- length(uDupENSP)
noExactMatch <- character()

for (i in 1:N) {
  percIDs <- IDs$percID[IDs$ENSP == uDupENSP[i]] # fetch the %ID values
  if (max(percIDs) < 99.999) {                   # if best is less than 100% ...
    noExactMatch <- c(noExactMatch, uDupENSP[i]) # ... store the ENSP ID
  }
}

# 36 proteins ... Let's look at the first:
IDs[IDs$ENSP == noExactMatch[1], ]
# ENSP00000403270 has UniProt IDs Q5T2P9 and A6NIR3. According to UniProt,
# Q5T2P9 is an obsolete entry, and A6NIR3 actually maps to a different STRING
# ENSP:
IDs[IDs$ENSP == "ENSP00000363207", ]

# Obsolete entries are removed from UniProt when gene-model errors have been
# corrected, or the genes were determined to be pseudogenes. We see that it
# would be incorrect to choose a similar UniProtID if an entry in UniProtKB is
# obsoleted, because that ID might actually be an exact match of a different
# ENSP.

# Therefore we should:
#   a: remove the 36 entries that have no exact match in UniProt,
#   b: keep only one (the first) entry of 100% ID for all duplicated ENSPs.

# a: removing the duplicated entries without exact match:
sel <- numeric()
for (i in 1:length(noExactMatch)) {
  sel <- c(sel, which(IDs$ENSP == noExactMatch[i]))
}
IDs <- IDs[-sel, ]   # 122 rows

# b: removing all but the first perfect match from duplicates

# First, update the list of unique, duplicate ENSP IDs, since we have removed
# some:
uDupENSP <- unique(IDs$ENSP[which(duplicated(IDs$ENSP))])

# Then select all the rows in which such duplicated IDs occur, except for the
# first one of each with 100% ID betwen the UniProt entry and the ENSP protein:
sel <- numeric()
for (i in 1:length(uDupENSP)) {
  myRows <- which(IDs$ENSP == uDupENSP[i])
  myRows <- myRows[-(which(IDs$percID[myRows] == 100)[1])]   # complicated :-)
  sel <- c(sel, myRows)
}

# ... and delete them from IDs
IDs <- IDs[-sel, ]   # remove 1,303 rows

# confirm:
which(duplicated(IDs$ENSP)) # none

# But ...
dupUni <- which(duplicated(IDs$UniProt))
length(dupUni) # 471

# Now, what about these 471 duplicated UniProt IDs?

IDs[IDs$UniProt == IDs$UniProt[which(duplicated(IDs$UniProt))[1]], ]
# Uniprot lists Q96P26 as NT5C1B, and crossreferences only ENSP00000352904.
# STRING db returns "NT5C1B	5’-nucleotidase, cytosolic IB (610 aa)" for
# ENSP00000352904, but "ENSG00000250741	NT5C1B-RDH14 readthrough (602 aa)" for
# ENSP00000433415 (note that the ENSP IDs are not identical). Both the latter
# ENSP IDs identify a fusion protein, which is described in the NCBI Gene
# database as a "...naturally occurring read-through transcription between the
# neighboring NT5C1B (5'-nucleotidase, cytosolic IB) and RDH14 (retinol
# dehydrogenase 14) genes on chromosome 2. Alternative splicing results in
# multiple transcript variants, one of which encodes a fusion protein that
# shares sequence identity with the products of each individual gene". The
# interaction edges between the two ENSP IDs in STRING have no overlap, simply
# creating a union of the two sets of edges under the common UniProt identifier
# would not be justified. I think that this is a case that cannot be
# automatically resolved.
#
tsv[grep("Q96P26", tsv$V2), ]

# In the case of ...
IDs[IDs$UniProt == IDs$UniProt[which(duplicated(IDs$UniProt))[2]], ]
# ... both STRING entries point to the same protein - once listed by its symbol,
# and once by its ENSG identifier. Both STRING entries have the same
# interactors. This is a case that also can't be resolved automatically, without
# considering additional data (e.g. whether the interaction partners are
# identical.)

# In the case of ...
IDs[IDs$UniProt == IDs$UniProt[which(duplicated(IDs$UniProt))[3]], ]
# ... two different ENSP identifiers point to the same UniProt ID - but the ENSP
# proteins don't have the same sequence. This case could be resolved by deleting
# rows with duplicated UniProt IDs that don't have 100% sequence identity.

# In the case of ...
IDs[IDs$UniProt == IDs$UniProt[which(duplicated(IDs$UniProt))[4]], ]
# ... the STRING entries point to two different proteins, with non-overlapping
# interactors, only the first of which corresponds to the UniProt KB entry for
# the UniProt ID. Again, this case could not be automatically resolved.

# etc.  You are welcome to explore additional examples by using the  UniProtKB,
# STRING and HGNC resources, to ask yourself how _you_ would resolve these
# issues.

# In summary, there are a number of reasons for duplicated UniProt IDs, many of
# which appear to point to uncertain annotations, if not outright error. I think
# it is therefore prudent to remove rows with duplicate UniProt IDs entirely.
#
# Note that it is not sufficient to remove the rows we have listed in dupUni,
# since duplicated()  returns only the second of a duplicated pair, not both.
# Consider:
which(duplicated(c(1, 2, 3, 1)))

# We need to compile an explicit selection ...

sel <- numeric()
for (i in 1:length(dupUni)) {
  sel <- c(sel, which(IDs$UniProt == IDs$UniProt[dupUni[i]]))
}
# ... and remove the offending rows.
IDs <- IDs[-sel, ]   # 1,458 rows

# Confirm:
which(duplicated(IDs$ENSP))
which(duplicated(IDs$UniProt))

# ... we now have a unique mapping of ENSP IDs to UniProt IDs - although it is
# not exhaustive. How many ENSP IDs did we have to give up?
( a <- length(unique(tsv$V3)) )  # 19,110
( b <- nrow(IDs)              )  # 18,216
(100 * (a - b)) / a              # 4.7 %


# ==============================================================================
#        PART TWO: UNIPROT TO HGNC SYMBOL
# ==============================================================================

# Next - we wanted to find gene symbols. We write a text-file that we can upload
# to the UniProt ID mapping service ...
writeLines(IDs$UniProt, "UniProtIDs.txt")
# ... and we can download and read the result: "18,144 out of 18,216 identifiers
# from UniProtKB AC/ID were successfully mapped to 18,150 Gene name IDs." I have placed the resulting file into the data directory (again: compressed).

uni2symMap <- read.delim("data/UniProt-SymbolMap.2017-02-01.txt.gz",
                  header = TRUE,
                  stringsAsFactors = FALSE)
colnames(uni2symMap) <- c("UniProt", "sym")
head(uni2symMap)
nrow(uni2symMap)  # 18,155

# Can we now use these symbols to translate UniProt IDs? Not before we check if
# they are unique ...
which(duplicated(uni2symMap$UniProt))
which(duplicated(uni2symMap$sym))
# Hm. Here we go again.

# Duplications in the $UniProt column mean: one UniProt ID is mapped to several
# symbols. We'll need to see how these are related.

dupUni <- which(duplicated(uni2symMap$UniProt))
for (i in 1:length(dupUni)) {
  iDup <- which(uni2symMap$UniProt == uni2symMap$UniProt[dupUni[i]])
  for (j in 1:length(iDup)) {
    cat(sprintf("%d\t%s\t%s\n",
                iDup[j],
                uni2symMap[iDup[j], "UniProt"],
                uni2symMap[iDup[j], "sym"]))
  }
  cat("\n")
}
# I think it is reasonable that these isoforms could be collapsed to the protein
# referenced in the first entry i.e. deleting all rows that are referenced in
# dupUni.
uni2symMap <- uni2symMap[-(dupUni), ]

# Duplications in the $sym column mean: different UniProt IDs are mapped to the
# same symbol.
dupSym <- which(duplicated(uni2symMap$sym))
for (i in 1:length(dupSym)) {
  iDup <- which(uni2symMap$sym == uni2symMap$sym[dupSym[i]])
  for (j in 1:length(iDup)) {
    cat(sprintf("%d\t%s\t%s\n",
                iDup[j],
                uni2symMap[iDup[j], 1],
                uni2symMap[iDup[j], 2]))
  }
  cat("\n")
}
# These Uniprot IDs need to be manually checked in UniProt, by retrieving their
# UniProt KB pages - eg: http://www.uniprot.org/uniprot/P0CG13 and
# http://www.uniprot.org/uniprot/P0CG13.

# No STRING interactors for P0CG12.
# No STRING interactors for Q9HC47.
# No STRING interactors for Q8N9K7.

# Annotating P86397 and O95059 to the symbol RPP14 is a curation error. These
# are two entirely different proteins (the former a member of the ribonuclease P
# complex, the latter a mitochondrial thioester dehydratase (HTD2)) - but while
# they are transcribed from the same genomic location, this is in fact a
# polycistronic transcript resulting in independent translation
# (https://www.ncbi.nlm.nih.gov/pubmed/17898086). Nevertheless, they are both
# assigned the RPP14 gene symbol. We will use HTD2 for P86397 instead.

# Q9H496 is an isoform of Q8NFQ8 (TOR1AIP2) that is annotated for only a subset
# of interactions of Q8NFQ8. We can delete it.

# These cases demonstrate how difficult it is to accommodate the large spectrum
# of biological facts in a database schema. We have unclear mapping from DNA to
# mRNA and from mRNA to protein, we need to account for database updates, and
# updates of genome assembly, we need to take isoforms and perhaps PTMs into
# account ... yet, we also need simple identifiers that allow us to
# cross-reference knowledge. ENSP IDs are too specific, UniProt IDs are too
# abstract, gene symbols are not always specific enough, still they appear to
# provide the best compromise. At least, the number of duplicated gene symbols
# is small enough that the last set of changes can be made manually.

sel <- which(uni2symMap$UniProt %in% c("P0CG12", "Q9HC47", "Q8N9K7", "Q9H496"))
uni2symMap <- uni2symMap[-sel, ]
sel <- which(uni2symMap$UniProt == "P86397")
uni2symMap$sym[sel] <- "HTD2"

# Confirm ...
which(duplicated(uni2symMap$UniProt))
which(duplicated(uni2symMap$sym))

# Note: while we have a usable correspondence of IDs, it is not perfect because
# not all ENSP IDs were mapped in the first place, but more importantly, some of
# what we did required manual intervention - and this means we have a workflow
# that won't scale. How to address such a situation is an open question.


# ==============================================================================
#        PART THREE: THE MAP
# ==============================================================================


# We have now a unique mapping of UniProt IDs to gene symbols. We just
# need to merge the two tables, using the UniProt IDs. We could write a for-loop
# for this, but there is also a specific R function: merge()
?merge

ID2symMap <- merge(IDs, uni2symMap,
                   by.x = "UniProt", by.y = "UniProt",
                   all.x = TRUE)

head(ID2symMap)

# Confirm:
i <- 1234
ID2symMap[i, ]
IDs[IDs$UniProt == ID2symMap$UniProt[i], ]
uni2symMap[uni2symMap$UniProt == ID2symMap$UniProt[i], ]

i <- 2345
ID2symMap[i, ]
IDs[IDs$UniProt == ID2symMap$UniProt[i], ]
uni2symMap[uni2symMap$UniProt == ID2symMap$UniProt[i], ]

# Finishing up: two things remain to be done. First of all, some Uniprot IDs
# could not be matched to gene symbols. We will remove these rows from the
# mapping table. We could have perhaps used UniProt IDs for these, but we would
# then be mixing the semantics of our identifiers and that may cause problems
# later on, in case we attempt to cross-reference information based on the
# identifier.

sel <- which(is.na(ID2symMap$sym))  # 76 rows
ID2symMap <- ID2symMap[-sel, ]

# Secondly, in order to use the map efficiently, we need a convenient way to
# search entries. This is best done by using the ENSP IDs as rownames. Therefore
# we don't actually need the ENSP column anymore. Nor do we need the percID
# column - and while we are at it, we might split the map into two, to map to
# symbols or UniProt IDs, This leaves us with only a single column each, thus we
# don't need this to be a data frame anymore but can convert it to a named
# vector:
tmp <- ID2symMap
ID2symMap <- tmp$sym
names(ID2symMap) <- tmp$ENSP
head(ID2symMap)

ID2UniMap <- tmp$UniProt
names(ID2UniMap) <- tmp$ENSP
head(ID2UniMap)

rm(tmp)

# How do we use these ? Simple:

for (i in 1:10) {
  ENSP <- sample(IDs$ENSP, 1)
  cat(sprintf("%s\t%10s\t%s\n",ENSP, ID2symMap[ENSP], ID2UniMap[ENSP]))
}

# Is it fast? Let's try 100,000 mappings.
tmp <- sample(IDs$ENSP, 100000, replace = TRUE)
head(tmp)

# Now map all of these to symbols.
tmp2 <- ID2symMap[tmp]  # ... pretty much instantaneous.
head(tmp2)

# Compare this with a for-loop, mapping only the first 10000 elements of tmp and
# searching for the matches in a dataframe:
slowMap <- data.frame(ID = names(ID2symMap),
                      sym = ID2symMap,
                      stringsAsFactors = FALSE)
N <- 10000
tmp2 <- character(N)
for (i in 1:N) {
  sym <- slowMap$sym[which(slowMap$ID == tmp[i])]
  if (length(sym) > 0) {
    tmp2[i] <- sym
  }
}

# Final cleanup:
save(ID2symMap, file = "ID2symMap.RData")
save(ID2UniMap, file = "ID2UniMap.RData")

rm(list = c("IDs", "uni2symMap","slowMap", "tsv",
            "a", "b", "i", "j", "myRows", "N", "sym",
            "dupENS", "dupSym", "dupUni", "uDupENSP", "iDup",
            "noExactMatch", "percIDs", "sel", "tmp", "tmp2"))

# To use the mapping files elsewhere, simply load() them from this directory.
# Remember that you will need to use an explicit path to this project directory
# to be able to load the files:
getwd()

# That's all.


# [END]
