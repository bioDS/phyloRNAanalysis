#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
# unofficial bash strict mode, see:
# http://redsymbol.net/articles/unofficial-bash-strict-mode/

# get_gatk_resources.sh
#
# downloads the following gatk resources into the current directory:
# pon -- 1000g_pon.hg38.vcf.gz + tbi
# gnomad -- af-only-gnomad.hg38.vcf.gz + tbi
gs="gs://gatk-best-practices/somatic-hg38"


gsutil cp "$gs/1000g_pon.hg38.vcf*" .
gsutil cp "$gs/af-only-gnomad.hg38.vcf.gz*" .
