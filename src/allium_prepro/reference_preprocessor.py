import os
import requests
import pandas as pd
import rpy2.robjects as robjects


class ReferencePreprocessor():
    def __init__(self,
                 genome_version='Homo_sapiens.GRCh38.103'):
        self._genome_version = genome_version

        # Ensure the ref directory exists
        self._ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     '../data/reference')
        os.makedirs(self._ref_dir, exist_ok=True)

        # Set filenames
        self._annot_file_full = \
            os.path.join(self._ref_dir,
                         f'{genome_version}.allium.annotations.full.csv')
        self._annot_file_filtered = \
            os.path.join(self._ref_dir,
                         f'{genome_version}.allium.annotations.filtered.csv')

    def run(self):
        self._download_ref()
        self._parse_gtf()
        self._filter()
        self._cleanup()

    def _download_ref(self):
        REFERENCE_GENOME = f'{self._genome_version}.gtf.gz'
        REFERENCE_GENOME_URL = f'http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/{REFERENCE_GENOME}'

        # If reference file already exists, say so and exit
        self._gtf_file_path = os.path.join(self._ref_dir, REFERENCE_GENOME)

        if os.path.exists(self._gtf_file_path):
            print(f"{REFERENCE_GENOME} already exists at {self._gtf_file_path}. Done.")
            return

        print(f"Downloading reference from {REFERENCE_GENOME_URL} to {self._ref_dir}...")
        # Download the file
        response = requests.get(REFERENCE_GENOME_URL, stream=True)
        response.raise_for_status()  # Check if the request was successful

        # Save the file to the specified directory
        with open(self._gtf_file_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        print(f"Downloaded {REFERENCE_GENOME}.gz to {self._gtf_file_path}")

    def _parse_gtf(self):
        print('Parsing GTF...')

        robjects.r('''
            library(Rsubread)
            library(sva)
            library(ballgown)

            parse_gtf <- function(gtf_file_path, output_file_path) {
                SAF <- Rsubread::flattenGTF(gtf_file_path)
                GeneLength <- rowsum(SAF$End-SAF$Start+1, SAF$GeneID)
                annot <- gffRead(gtf_file_path)

                # Retrieve the attributes
                annot$gene_id = getAttributeField(annot$attributes,
                    field = "gene_id")

                annot$gene_name = getAttributeField(annot$attributes,
                    field = "gene_name")

                annot$gene_biotype = getAttributeField(annot$attributes,
                    field = "gene_biotype")

                tokeep <- c("seqname", "feature", "gene_id", "gene_name", "gene_biotype")
                annot_filtered <- annot[tokeep]

                # Remove quotes
                annot_filtered[] <- lapply(annot_filtered, gsub, pattern='"', replacement='')

                # Keep only genes
                annot_filtered <- annot_filtered[annot_filtered$feature == 'gene', ]

                # Convert gene lengths to dataframe
                gene_lengths <- as.data.frame(GeneLength)
                gene_lengths$gene_id <- rownames(GeneLength)

                # Reorder columns to have gene_id as the first column
                gene_lengths <- gene_lengths[, c("gene_id", names(gene_lengths)[-ncol(gene_lengths)])]

                # Name the columns
                colnames(gene_lengths) <- c("gene_id", "gene_length")

                # Join DataFrames on 'gene_id' column
                full_annot <- merge(annot_filtered, gene_lengths, by = "gene_id", all = TRUE)

                # Rename columns
                colnames(full_annot) <- c("id", "chr", "feature", "name", "biotype", "length")

                # Drop feature column
                full_annot$feature <- NULL

                # Strip trailing semicolons from the biotype column
                full_annot$biotype <- sub(";$", "", full_annot$biotype)

                # Sort annotations by id
                full_annot <- full_annot[order(full_annot$id), ]

                # Write processed metadata
                write.table(full_annot, output_file_path, col.names = TRUE, row.names = FALSE, sep = ',', quote = FALSE)
            }
        ''')

        r_parse_gtf = robjects.r['parse_gtf']
        r_parse_gtf(self._gtf_file_path, self._annot_file_full)
        print(f"Created {self._annot_file_full}")

    def _filter(self):
        # Read in file
        annot = pd.read_csv(self._annot_file_full)

        # Get all rows where biotype is protein_coding
        protein_coding_annot = annot.loc[annot['biotype'] == 'protein_coding']

        # Only keep genes that are on chrs 1-22 and X
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['X']

        # Keep only genes where chr is in valid_chromosomes
        filtered_annot = protein_coding_annot[
            protein_coding_annot['chr'].isin(valid_chromosomes)]

        # Remove ribosomal genes
        filtered_annot = filtered_annot[
            ~filtered_annot['name'].str.startswith(('RPL', 'RPS'))]

        # Write the filtered DataFrame to a CSV file
        filtered_annot.to_csv(self._annot_file_filtered, index=False)

    def _cleanup(self):
        print("Cleaning up...")
        # Remove the original GTF file
        os.remove(self._gtf_file_path)
        print("Done.")
