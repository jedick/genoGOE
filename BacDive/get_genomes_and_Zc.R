# BacDive/get_genomes.R
# Get protein sequences for genome accessions
# 20260313 jmd

# Check for genomes directory
if(!dir.exists("genomes")) stop("Please create 'genomes' directory")

# Read BacDive data
dat <- read.csv("cleaned_data.csv", check.names = FALSE)
# Get GCA accessions
GCA_accessions <- dat$`sequence_genomes.INSDC accession`
# Make GCF accessions
GCF_accessions <- paste0(gsub("GCA", "GCF", GCA_accessions), ".1")

# Loop over GCF accessions
for(GCF in GCF_accessions) {
  # Skip missing accessions
  if(GCF == ".1") next
  # Check for existing zip file
  outfile <- file.path("genomes", paste0(GCF, ".zip"))
  if(!file.exists(outfile)) {
    print(outfile)
    # Download protein sequences using NCBI's datasets command
    cmd <- paste("datasets download genome accession", GCF, "--include protein")
    system(cmd)
    file.rename("ncbi_dataset.zip", outfile)
  }
}

# Create new data frame with empty Zc column
out <- dat
out$Zc <- NA
# Loop over rows
for(i in 1:nrow(out)) {
  # Look for file with name of GCF accession
  zipfile <- paste0("genomes/", GCF_accessions[i], ".zip")
  if(file.exists(zipfile)) {
    # Unzip the archive, omitting README.md and md5sum.txt
    cmd <- paste("unzip", zipfile, "-x README.md md5sum.txt")
    system(cmd)
    faafile <- file.path("ncbi_dataset/data/", GCF_accessions[i], "protein.faa")
    if(file.exists(faafile)) {
      aa <- canprot::read_fasta(faafile)
      sumaa <- canprot::sum_aa(aa)
      # Calculate Zc and insert into output data frame
      Zc <- canprot::Zc(sumaa)
      out$Zc[i] <- round(Zc, 6)
    }
    # Remove unzipped directory
    system("rm -rf ncbi_dataset/ README.md")
  }
}
# Save output
write.csv(out, "cleaned_data_with_Zc.csv", row.names = FALSE)


