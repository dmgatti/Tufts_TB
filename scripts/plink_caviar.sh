#!/vin/bash

PLINK_EXE='/home/dmgatti/software/plink_linux_x86_64_20201019/plink'
CAVIAR_EXE='/home/dmgatti/software/caviar/CAVIAR-C++/CAVIAR'

INPUT_DIR='/media/dmgatti/hdb/projects/TB/manuscripts/manuscript1/results/caviar'
OUTPUT_DIR='/media/dmgatti/hdb/projects/TB/manuscripts/manuscript1/results/caviar'

# Convert PED & MAP filed to BED/BIM/FAM files.
${PLINK_EXE} --file ${INPUT_DIR}/pc1 --mouse --make-bed --out ${OUTPUT_DIR}/pc1

# Get LD using r.
${PLINK_EXE} --bfile ${INPUT_DIR}/pc1 --mouse --r square --out ${OUTPUT_DIR}/pc1

# Get association mapping without kinship.
${PLINK_EXE} --bfile ${INPUT_DIR}/pc1 --mouse --linear standard-beta --covar ${OUTPUT_DIR}/pc1.covar --out ${OUTPUT_DIR}/pc1\
  --parameters 1

# Get Z-scores from PLINK file. (Not 100% sure these are Z-scores. Can't figure out how to get SE)
awk -F ' ' '{print $2, $7}' ${OUTPUT_DIR}/pc1.assoc.linear > ${OUTPUT_DIR}/pc1.zscore

# Run CAVIAR.
${CAVIAR_EXE} -l ${OUTPUT_DIR}/pc1.ld -z ${OUTPUT_DIR}/pc1.zscore -c 10 -o ${OUTPUT_DIR}/pc1.caviar
