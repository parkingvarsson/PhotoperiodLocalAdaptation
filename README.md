# PhotoperiodLocalAdaptation
Scripts for Wang et al. (2017) "Natural variation at FLOWERING LOCUS T2 mediates local adaptation in a key life history trait in European aspen". https://www.biorxiv.org/content/early/2017/12/13/178921

<h2>Documentation of Scripts</h2>

<b>1. Sequencing quality checking, read mapping and post-mapping filtering</b>

    QualityControl.sh - Use Trimmomatic v0.30 and FastQC to do sequence quality checking

    runBWA_mem_asp201.SwAsp.sh - Use BWA-MEM to do read mapping

    Picard_markduplicates_SwAsp.sh - Use Picard to correct for artifacts of PCR duplication

    realign.sh - Use GATK to do read realignment around indels

<b>2. SNP and genotype calling</b>

        HaplotypeCaller_mem_tremula.sh - Use GATK to do SNP calling

        snpEff.gatk.sh - Use snpEff to annotate the SNPs

        vcf_to_maf.sh vcf2maf.pl - Use perl script to map each variant to only one of all possible gene isoforms

<b>3. Relatedness, population structure and isolation-by-distance</b>

        plink_prunedLD.sh - Use PLINK to generate Linkage-disequilibrium(LD)-trimmed SNP sets

        eigen_pop_genetics.sh - Use smartpca program in EIGENSOFT to perform PCA analysis

        fst_matrix.R - create the matrix of Fst estimates

        Isolation_by_distance.R - create the matrix of geographic distance

        SwAsp.IBD.plot.R - plot Isolation-by-distance

        mantel.R - Mantel test 

<b>4. Screening for SNPs associated with local adaptation</b>

        PCAdapt.SwAsp94.R - PCAadapt test

        pcadapt_fst_cor.R - The relationship between PCadapt results and Fst values

        pca_corr_env.R - PCA analysis for the environment data

        env_PC.error_bar.R - Relationship between the environmental PC1 scores and the number of days with degree higher than 5

        LEA.K1.R LEA.K2.R LEA.K3.R - use a latent factor mixed-effect model (LFMM) implemented in the R package LEA to detect SNPs associated with first environmental PC, with the latent factors (K) from 1 to 3.

        LEA_zscore.94samples.R - Transform the z-scores from LFMM results to p-values

        Ekebo_Savar_budset.as - Use Asreml to estimate the genetic values of bud set

        GEMMA.lmm.noPC.budset.sh - Use GEMMA to do GWAS for budset

<b>5. Genotype imputation</b>

        define_ancestral.SwAsp.pseudo_chr.sh - Use BEAGLE to do genotype imputation and also define ancestal and derived allele based on the sequences of outgroup species of P.tremuloides and P.trichocarpa

        Beagle_comparison.sh - Create several missing genotypes for simulation

        Create_simu_missing_file.sh - Tor simulate by BEAGLE and calculate accuracy

        createFile_imputate.R - Create imputation files with different level of missing values

        impute_accuracy_beagle_vcf.pl - Calculate the imputation accuracy by each sample and SNP

        impute_evaluate.R - Evaluate the imputation accuracy


<b>6. Positive selection</b>

        angsd_SFS_SwAsp.all.sh - Use ANGSD to estimate the genetic diversity in specific groups of populations

        genome_wide.summary.compare.sh - Use vcftools, selscan and H12 test to perform a set of selection tests across the genome

        select.chr10.summary.sh - Use vcftools, selscan and H12 test to perform a set of selection tests on the specific region (~700kbp) on Chr10

        sweepfinder2.FT2paper.sh - Run SweepFinder2 to detect selection signals

        ehh.plot.sh colormap.plotting.R - Create EHH plot

        plink.ldheatmap.sh LDheatmap.R - Use PLINK to calculate LD across SNPs

        chr10_genome.sig.angsd_tP_tajD.group_plot.R - Compare the genetic diversity between chr10 region and genome-wide levels

        chr10_genome.sig.fst.group_plot.R - Compare Fst between chr10 region and genome-wide levels

        chr10_genome.sig.h12.group_plot.R - Compare H12 and H2/H1 between chr10 region and genome-wide levels

        chr10_genome.sig.sweepfinder2.group_plot.R - Compare CLR between chr10 region and genome-wide levels

        chr10_genome.sig.ihs_nsl.R - Compare iHS and nSL values between chr10 region and genome-wide levels

        caviar.chr10.sh caviar.z-score.R - Run CAVIAR on specific region of Chr10

        run_sweep_sim.sg - Simulate independent selective sweep events using the coalescent simulation program msms and analyse the results using SweepFinder2 to assess region affected by sweep

        ms2sf2.pl - Convert msms output to input format for SweepFinder2

        ABCinference.R - Script for performing Approximate Bayesian Computation (ABC) to jointly estimate s (the strength of selection on the beneficial mutation causing the sweep) and T (the time since the beneficial allele fixed). Script modified from original provided by Ormond et al (2016) at http://jjensenlab.org/wp-content/uploads/2016/02/ABC_inference.zip

