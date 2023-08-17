#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time 72:00:00
#SBATCH --job-name seed
#SBATCH -e log.seed.err
#SBATCH -o log.seed.out

module load R/3.6.1

covars_label="optional"
model=FNW1_svmRadial
scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

if [ ! -f ${outdir}rds/${id}lfs_TSS200_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate TSS200 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate TSS200 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_TSS1500_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate TSS1500 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate TSS1500 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_3UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate 3UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate 3UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_5UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate 5UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate 5UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_1stExon_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R  --outdir $outdir --id $id --lfs --aggregate 1stExon $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate 1stExon $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_gene_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --gene $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --gene $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}lfs_Body_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate Body $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate Body $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_TSS200_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate TSS200 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate TSS200 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_TSS1500_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate TSS1500 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate TSS1500 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_1stExon_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate 1stExon $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate 1stExon $covars $pred_covars
fi             
if [ ! -f ${outdir}rds/${id}nocancer_gene_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --gene $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --gene $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_Body_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate Body $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate Body $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_3UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate 3UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate 3UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_5UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --cancer --aggregate 5UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --cancer --aggregate 5UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS200_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate TSS200 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate TSS200 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS1500_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate TSS1500 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate TSS1500 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_1stExon_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate 1stExon $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate 1stExon $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_gene_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --gene $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --gene $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_Body_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate Body $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate Body $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_3UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate 3UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate 3UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}nocancer_lfs_5UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --cancer --aggregate 5UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --cancer --aggregate 5UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}gene_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --gene --outdir $outdir --id $id $covars
	Rscript $scripts_dir/predictSingleExtVal.R --gene --outdir $outdir --id $id $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}Body_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --aggregate Body --outdir $outdir --id $id $covars
	Rscript $scripts_dir/predictSingleExtVal.R --aggregate Body --outdir $outdir --id $id $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}TSS200_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --aggregate TSS200 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --aggregate TSS200 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}TSS1500_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --aggregate TSS1500 $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --aggregate TSS1500 $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}3UTR_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --aggregate 3UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --aggregate 3UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}5UTR_${covars_label}${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --aggregate 5UTR $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --aggregate 5UTR $covars $pred_covars
fi
if [ ! -f ${outdir}rds/${id}1stExon_${covars_label}_${model}_ageofonset_ROCInfoTest.rds ] ; then
	Rscript $scripts_dir/runSingleModel.R  --outdir $outdir --id $id --aggregate 1stExon $covars
	Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --aggregate 1stExon $covars $pred_covars
fi

