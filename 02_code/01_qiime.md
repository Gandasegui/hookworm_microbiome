# Code to import the Fastq files, read and generate QZAs

### Generate work environment

```bash
mkdir /home/jgandasegui/R/MARS_hookworm/NEW_IMPORT/
cd /home/jgandasegui/R/MARS_hookworm/NEW_IMPORT/
mkdir qiime_art
```

### Generate manifest.tsv doc for importaitn the secuences
```R
path<-"~/R/MARS_hookworm/NEW_IMPORT"
run<-"fastq_files"

forward<-list.files(glue("{path}/{run}"), pattern = "R1_001.fastq.gz", full.names = TRUE)
reverse<-list.files(glue("{path}/{run}"), pattern = "R2_001.fastq.gz", full.names = TRUE)

##Comprobamos que las muestras forward y reverse estén en el mismo orden
all.equal(str_replace(reverse,  "_L001_R2_001.fastq.gz", ""),
          str_replace(forward, "_L001_R1_001.fastq.gz", ""))

manifest<-data.frame(str_replace(reverse, "_L001_R2_001.fastq.gz", ""))%>%rename(`sample-id`=1)%>%
  as_tibble()%>%
  mutate(`sample-id`=str_replace(`sample-id`,paste0(path,"/", run, "/Fastq/"),""),
         `forward-absolute-filepath`=forward,
         `reverse-absolute-filepath`=reverse)%>%
  filter(`sample-id`!="Undetermined_S0")

manifest$`sample-id` <- sub("_[^_]+$", "", manifest$`sample-id`)
manifest$`sample-id` <- gsub("^[^:]*/", "", manifest$`sample-id`)

write.table(manifest, file = "manifest.tsv", row.names=FALSE, sep="\t", quote = FALSE)
```

### Now qiime2 is used
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path qiime_art/demux-paired-end_hook.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data qiime_art/demux-paired-end_hook.qza \
  --o-visualization qiime_art/demux-paired-end_hook.qzv
qiime tools view qiime_art/demux-paired-end_hook.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime_art/demux-paired-end_hook.qza \
  --p-trim-left-f 30 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 220 \
  --o-representative-sequences qiime_art/dada_rep-seqs.qza \
  --o-table qiime_art/dada_table.qza \
  --o-denoising-stats qiime_art/dada_stats.qza

# To see whatever...
qiime metadata tabulate \
  --m-input-file qiime_art/dada_stats.qza \
  --o-visualization qiime_art/dada_stats.qzv
qiime tools view qiime_art/dada_stats.qzv

# Summarises fueatures table samples
qiime feature-table summarize \
  --i-table qiime_art/dada_table.qza \
  --o-visualization qiime_art/dada_table.qzv \
  --m-sample-metadata-file metadata_mod.tsv #firlt coluimn sample-id and conver to tsv
qiime tools view qiime_art/dada_table.qzv

# Summarises 
qiime feature-table tabulate-seqs \
  --i-data qiime_art/dada_rep-seqs.qza \
  --o-visualization qiime_art/dada_rep-seqs.qzv
qiime tools view qiime_art/dada_rep-seqs.qzv

##Let's generate the rooted tree - phygenetics
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime_art/dada_rep-seqs.qza \
  --o-alignment qiime_art/aligned-rep-seqs.qza \
  --o-masked-alignment qiime_art/masked-aligned-rep-seqs.qza \
  --o-tree qiime_art/unrooted-tree.qza \
  --o-rooted-tree qiime_art/rooted-tree.qza

#Now, repeat with the silva classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads qiime_art/dada_rep-seqs.qza \
  --o-classification qiime_art/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiime_art/taxonomy.qza \
  --o-visualization qiime_art/taxonomy.qzv
qiime tools view qiime_art/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime_art/dada_table.qza \
  --i-taxonomy qiime_art/taxonomy.qza \
  --m-metadata-file metadata_mod.tsv \
  --o-visualization qiime_art/taxa-bar-plots_silva.qzv
qiime tools view qiime_art/taxa-bar-plots_silva.qzv

mv dada_table.qza table.qza
```
### Now I load into R:
-rooted-tree.qza

-metadata_mod.tsv

-table.qza

-taxonomy.qza
