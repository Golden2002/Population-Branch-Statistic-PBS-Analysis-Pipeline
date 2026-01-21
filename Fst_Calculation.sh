#!/bin/bash
#SBATCH --job-name=Population_Analysis
#SBATCH --output=/path/to/logs/Analysis_%j.log
#SBATCH --error=/path/to/logs/Analysis_err_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --nodes=1

######################################
# Initialize conda environment
######################################
__conda_setup="$('/path/to/miniconda3/bin/conda' 'shell.bash' 'hook' 2>/dev/null)"
if [ $? -eq 0 ]; then eval "$__conda_setup"
else
  if [ -f "/path/to/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/path/to/miniconda3/etc/profile.d/conda.sh"
  else
    export PATH="/path/to/miniconda3/bin:$PATH"
  fi
fi
unset __conda_setup
conda activate genome || { echo "❌ Conda environment activation failed"; exit 1; }

######################################
# Parameter Settings
######################################
VCF="/path/to/data/chr_Auto.vcf.gz"
BASE_DIR="/path/to/analysis_directory"
POP_LIST="${BASE_DIR}/population_list.txt"
POP_LABELS="${BASE_DIR}/population_labels.txt"
POP_DIR="${BASE_DIR}/population_samples"
WINDOW_SIZE=10000
WINDOW_STEP=5000
TAJIMA_WINDOW=5000   # Smaller window for higher resolution
TARGET_POP="TargetPopulation"

FST_DIR="${BASE_DIR}/fst_sliding_window"
TARGET_FST_DIR="${BASE_DIR}/target_fst"
TAJIMA_DIR="${BASE_DIR}/tajimaD"
MATRIX_OUT="${BASE_DIR}/fst_matrix.txt"

mkdir -p $FST_DIR $TARGET_FST_DIR $TAJIMA_DIR $POP_DIR || exit 1

######################################
# Index VCF
######################################
if [ ! -f "${VCF}.tbi" ]; then
  echo "Creating index for VCF file..."
  bcftools index -t "$VCF" || { echo "❌ Indexing failed"; exit 1; }
fi

######################################
# Create sample lists for each population based on population labels
######################################
awk '{print $1 > "'${POP_DIR}'/"$2"_1col.txt"}' $POP_LABELS

######################################
# 1. Sliding Window Fst Parallel Calculation
######################################
module load parallel/20220922 vcftools/0.1.16

combs=$(mktemp)
awk 'NR==FNR{a[NR]=$1;next}{for(i=1;i<FNR;i++)print a[i],$1}' $POP_LIST $POP_LIST > $combs

fst_calc() {
  pop1="$1"
  pop2="$2"
  outfile="${FST_DIR}/${pop1}_vs_${pop2}"
  if [[ -s ${outfile}.windowed.weir.fst ]]; then
    echo "✔️ ${pop1} vs ${pop2} already exists, skipping"
  else
    vcftools --gzvcf "$VCF" \
             --weir-fst-pop "${POP_DIR}/${pop1}_1col.txt" \
             --weir-fst-pop "${POP_DIR}/${pop2}_1col.txt" \
             --fst-window-size $WINDOW_SIZE \
             --fst-window-step $WINDOW_STEP \
             --out "$outfile"
  fi
}
export -f fst_calc
export VCF POP_DIR FST_DIR WINDOW_SIZE WINDOW_STEP

echo "===== Parallel calculation of all Fst pairs ====="
parallel --colsep ' ' -j 16 fst_calc {1} {2} :::: $combs
rm $combs
echo "✅ Sliding window Fst calculation completed"

######################################
# 2. Extract Target Population Fst
######################################
# Uncomment and modify as needed
# for ref in "${REF_POPS[@]}"; do
#   src="${FST_DIR}/${TARGET_POP}_vs_${ref}.windowed.weir.fst"
#   dest="${TARGET_FST_DIR}/${TARGET_POP}_vs_${ref}.windowed.weir.fst"
#   if [ -f "$src" ]; then
#     cp "$src" "$dest"
#   else
#     echo "⚠️ ${src} missing, skipping"
#   fi
# done

######################################
# 3. Tajima D Analysis
######################################
TAJIMA_OUT="${TAJIMA_DIR}/${TARGET_POP}_TajimaD"
if [ -f "${TAJIMA_OUT}.Tajima.D" ]; then
  echo "✔️ Tajima D already exists"
else
  vcftools --gzvcf "$VCF" \
           --keep "${POP_DIR}/${TARGET_POP}_1col.txt" \
           --TajimaD $TAJIMA_WINDOW \
           --out "$TAJIMA_OUT"
fi

######################################
# 4. Output Average Fst Matrix (can be used for NJ tree)
######################################
# Uncomment to use
# echo "Outputting average Fst matrix between populations..."
# {
#   echo -ne "\t"
#   paste -sd'\t' $POP_LIST
#   while read p1; do
#     echo -ne "$p1"
#     while read p2; do
#       file="${FST_DIR}/${p1}_vs_${p2}.windowed.weir.fst"
#       if [ -f "$file" ]; then
#         mean_fst=$(awk '$4 != "NaN" {sum += $4; n++} END{if(n>0) print sum/n; else print "NA"}' "$file")
#       else
#         mean_fst="NA"
#       fi
#       echo -ne "\t$mean_fst"
#     done < $POP_LIST
#     echo
#   done < $POP_LIST
# } > "$MATRIX_OUT"

######################################
# Completion Notification
######################################
echo -e "\n✅ Analysis completed! Result directories:"
echo "Fst window results: $FST_DIR"
echo "Target population Fst extraction: $TARGET_FST_DIR"
echo "Tajima D results: $TAJIMA_DIR"
echo "Average Fst matrix: $MATRIX_OUT"
