#!/bin/bash

# conda env
#conda activate esviritu_postproc;

# Define variables
SAMPLE_LIST=$1 #pXXXX_samplesheet.tsv used for EsViritu

# create sample list
cut -f1 ${SAMPLE_LIST} | tail -n+2  > samples.txt

SAMPLES_FILE="samples.txt"

process_contig() {
  SAMPLE_ID=$1
  CONTIG_ID=$2
  CONTIG_SEQ=$3
  
  echo "Processing contig: ${CONTIG_ID}"
  
  THREADS=8
  # Use contig ID for coverage file lookup
  COVERAGE_INFO="${SAMPLE_ID}.virus_coverage_windows.tsv"
  OUTDIR="./assembly_stats/${SAMPLE_ID}/${CONTIG_ID}"

  # Make outdir
  mkdir -p ${OUTDIR}

  # Create a temporary fasta file for this contig
  CONTIG_FASTA="${OUTDIR}/${CONTIG_ID}_contig.fasta"
  echo ">${CONTIG_ID}" > ${CONTIG_FASTA}
  echo "${CONTIG_SEQ}" >> ${CONTIG_FASTA}

  ASSEMBLY=${CONTIG_FASTA}

  # Quality stats
  echo "Generating quality statistics..."

  STATS_FILE="${OUTDIR}/${CONTIG_ID}_detailed_stats_EVassembly.txt"

  # Create summary report
  echo "================================================" > ${STATS_FILE}
  echo "Assembly Quality Report: ${CONTIG_ID}" >> ${STATS_FILE}
  echo "================================================" >> ${STATS_FILE}
  echo "" >> ${STATS_FILE}

  # Assembly metrics
  echo "ASSEMBLY METRICS" >> ${STATS_FILE}
  echo "----------------" >> ${STATS_FILE}

  # Get assembly length and N count
  ASSEMBLY_LENGTH=$(grep -v ">" ${ASSEMBLY} | tr -d '\n' | wc -c)
  N_COUNT=$(grep -v ">" ${ASSEMBLY} | tr -d '\n' | grep -o 'N' | wc -l)
  PERCENT_N=$(awk "BEGIN {printf \"%.2f\", ($N_COUNT/$ASSEMBLY_LENGTH)*100}")
  COMPLETENESS=$(awk "BEGIN {printf \"%.2f\", (($ASSEMBLY_LENGTH-$N_COUNT)/$ASSEMBLY_LENGTH)*100}")

  echo "Assembly Length: ${ASSEMBLY_LENGTH} bp" >> ${STATS_FILE}
  echo "Number of Ns: ${N_COUNT}" >> ${STATS_FILE}
  echo "Percent Ns: ${PERCENT_N}%" >> ${STATS_FILE}
  echo "" >> ${STATS_FILE}

  # Coverage metrics
  echo "COVERAGE METRICS" >> ${STATS_FILE}
  echo "----------------" >> ${STATS_FILE}

  # Calculate average coverage from depth file (if exists)
  if [ -f "${COVERAGE_INFO}" ]; then
    # Extract just the accession part (remove sample prefix)
    # If CONTIG_ID is "p1389-EH1RI0246_KU950540.1", this extracts "KU950540.1"
    ACCESSION=$(echo "${CONTIG_ID}" | sed "s/^${SAMPLE_ID}_//")

    # Filter coverage data to only this contig's accession
    CONTIG_COV_DATA=$(awk -v acc="${ACCESSION}" '$1 == acc {print $5}' "${COVERAGE_INFO}")
    
    if [ -n "$CONTIG_COV_DATA" ]; then
      AVG_COV=$(echo "$CONTIG_COV_DATA" | awk '{sum+=$1; count++} END {printf "%.2f", sum/count}')
      echo "Average Coverage: ${AVG_COV}x" >> ${STATS_FILE}

      # Additional coverage stats for THIS contig only
      echo "$CONTIG_COV_DATA" | sort -n | awk '
        {
          a[NR]=$1
          count++
        }
        END {
          min = a[1]
          max = a[count]
          median = (count%2==1) ? a[int(count/2)+1] : (a[count/2] + a[count/2+1])/2;
          printf "Minimum Coverage: %.0f\n", min;
          printf "Maximum Coverage: %.0f\n", max;
          printf "Median Coverage: %.2f\n", median
        }
      ' >> ${STATS_FILE}

      # Positions with coverage < 10x for THIS contig
      LOW_COV=$(echo "$CONTIG_COV_DATA" | awk '$1 < 10' | wc -l)
      CONTIG_POSITIONS=$(echo "$CONTIG_COV_DATA" | wc -l)
      PERCENT_LOW=$(awk "BEGIN {printf \"%.2f\", ($LOW_COV/$CONTIG_POSITIONS)*100}")
      echo "Positions with <10x coverage: ${LOW_COV} (${PERCENT_LOW}%)" >> ${STATS_FILE}
    else
      echo "No coverage data found for accession: ${ACCESSION}" >> ${STATS_FILE}
      AVG_COV="NA"
    fi
  else
    echo "Coverage file not found: ${COVERAGE_INFO}" >> ${STATS_FILE}
    AVG_COV="NA"
  fi
  echo "" >> ${STATS_FILE}

  # Seqkit stats
  echo "SEQUENCE STATISTICS (seqkit)" >> ${STATS_FILE}
  echo "----------------------------" >> ${STATS_FILE}
  seqkit stats ${ASSEMBLY} >> ${STATS_FILE}
  echo "" >> ${STATS_FILE}

  echo "Detailed quality statistics saved to: ${STATS_FILE}"

  # Create TSV summary
  TSV_FILE="./assembly_stats/${CONTIG_ID}_summary_EVassembly.tsv"

  echo -e "ContigID\tasm_length\tN_num\tN_percent\tCompleteness\tavg_cov" > ${TSV_FILE}
  echo -e "${CONTIG_ID}\t${ASSEMBLY_LENGTH}\t${N_COUNT}\t${PERCENT_N}\t${COMPLETENESS}\t${AVG_COV}" >> ${TSV_FILE}

  echo "TSV summary saved to: ${TSV_FILE}"

  echo "=== DONE with ${CONTIG_ID}! ==="
}

process_sample() {
  SAMPLE_ID=$1
  
  echo "Processing sample: $SAMPLE_ID"
  SAMPLE=$SAMPLE_ID
  OUTDIR="./assembly_stats/${SAMPLE}"

  # Make outdir
  mkdir -p ${OUTDIR}
  
  # Find the assembly file
  ASSEMBLY="./${SAMPLE_ID}_final_consensus.fasta"
  
  if [ ! -f "$ASSEMBLY" ]; then
      echo "ERROR: Assembly file not found: $ASSEMBLY"
      # Write to no assembly file
      return 1
  fi

# Check if any header already has the sample ID prefix
if grep -q "^>${SAMPLE_ID}_" "${ASSEMBLY}"; then
    echo "Sample ID already in headers, skipping modification"
else
    echo "Adding sample ID to headers"
    sed -i "s/^>\(.*\)_consensus$/>${SAMPLE_ID}_\1/" "${ASSEMBLY}"
fi

  # Filter short contigs (<1000 bp non-N sequence)
  echo "Filtering contigs..."
  seqkit fx2tab "$ASSEMBLY" | \
  awk -F'\t' '{
    seq = $2;
    gsub(/[Nn]/, "", seq);
    if (length(seq) >= 1000) {
      print $0
    }
  }' | \
  seqkit tab2fx > "${OUTDIR}/${SAMPLE}_filtered.fasta"

  FILTERED_ASSEMBLY="${OUTDIR}/${SAMPLE}_filtered.fasta"
  
  # Check if any contigs passed filtering
  CONTIG_COUNT=$(grep -c "^>" "${FILTERED_ASSEMBLY}" 2>/dev/null || echo "0")
  
  if [ "$CONTIG_COUNT" -eq 0 ]; then
    echo "WARNING: No contigs passed filtering for ${SAMPLE_ID}"
    return 1
  fi
  
  echo "Found ${CONTIG_COUNT} contig(s) after filtering"

  # Parse the filtered fasta file and process each contig
  awk '/^>/ {
    if (seqlen) {
      print header"\t"seq
    }
    header=$0
    seq=""
    seqlen=0
    next
  }
  {
    seq=seq$0
    seqlen+=length($0)
  }
  END {
    if (seqlen) {
      print header"\t"seq
    }
  }' "${FILTERED_ASSEMBLY}" | while IFS=$'\t' read -r header sequence; do
    # Remove the ">" from header
    CONTIG_ID=$(echo "$header" | sed 's/^>//')
    
    # Process this contig
    process_contig "${SAMPLE_ID}" "${CONTIG_ID}" "${sequence}"
  done
}

# Export functions
export -f process_sample
export -f process_contig

# Run 3 jobs in parallel
parallel -j 3 --env process_sample --env process_contig process_sample :::: "$SAMPLES_FILE"

# Gzip all fasta files created during processing
echo ""
echo "Compressing all FASTA files..."
find ./assembly_stats -name "*.fasta" -type f -exec gzip {} \;
echo "FASTA compression complete!"

# Concatenate all TSV summary files
echo "Concatenating all summary statistics..."
COMBINED_TSV="./assembly_stats/All_contigs_summary.tsv"

# Get the first TSV file to extract header
FIRST_TSV=$(find ./assembly_stats -name "*_summary_EVassembly.tsv" | head -n 1)

if [ -n "$FIRST_TSV" ]; then
  # Copy header from first file
  head -n 1 "$FIRST_TSV" > "$COMBINED_TSV"
  
  # Append all data (skip headers) from all TSV files
  find ./assembly_stats -name "*_summary_EVassembly.tsv" -exec tail -n +2 {} \; >> "$COMBINED_TSV"
  
  echo "Combined summary saved to: $COMBINED_TSV"
  echo "Total contigs: $(tail -n +2 "$COMBINED_TSV" | wc -l)"
else
  echo "WARNING: No summary TSV files found!"
fi

# Find samples without assemblies
echo ""
echo "Checking for samples without assembly files..."
NO_ASSEMBLY_FILE="./assembly_stats/Samples_without_assembly.txt"
> "$NO_ASSEMBLY_FILE"  # Clear file

# In the end loop - uses lowercase local variable
while read -r sample_id; do  # Reading into a local variable
  if [ ! -f "./${sample_id}_final_consensus.fasta" ]; then  # Uses local variable
    echo "$sample_id" >> "$NO_ASSEMBLY_FILE"
  fi
done < "$SAMPLES_FILE"

# Report results
if [ -s "$NO_ASSEMBLY_FILE" ]; then
  NUM_NO_ASSEMBLY=$(wc -l < "$NO_ASSEMBLY_FILE")
  echo "========================================="
  echo "WARNING: $NUM_NO_ASSEMBLY sample(s) had no assembly file:"
  cat "$NO_ASSEMBLY_FILE"
  echo "See: $NO_ASSEMBLY_FILE"
  echo "========================================="
else
  echo "All samples had assembly files!"
  rm -f "$NO_ASSEMBLY_FILE"
fi