# Pipeline ITH

## Description

A pipeline to study intratumor heterogeneity (ITH) with Canopy<sup>[1]</sup>.

#### Marathon

This pipeline has been inspired from the Marathon pipeline<sup>[2]</sup> proposed by the authors of Canopy. Marathon is a description of a conceptual pipeline using Falcon and Canopy. It is not a functional and automated pipeline. The goal of the ITH pipeline is to propose a working and automated pipeline.

## General overview

<img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/pipeline_overview.png" style="width:70%; display:block; margin:auto;" />

## Notes

In this documentation, patient ID has been replaced with _##_, and tumors ID has been replaced with _T1_ and _T2_.

## Steps

#### Post-alignment

* tool : BWA
* input : a BAM file
* output : a BAM file
* scripts : scripts_cobalt/template_postalt.sh, scripts_cobalt/launch_postalt.sh

#### Germline calling

* tool : Platypus
* inputs : a normal BAM file, human genome reference file, regions file
* output : a normal VCF file
* scripts : scripts_cobalt/template_platypus.sh, scripts_cobalt/launch_platypus.sh

Then th VCF output file has been filtered on _PASS_ value : keep_pass.sh.

#### Somatic calling

* tool : Strelka2
* inputs : a tumor BAM file, its associated normal BAM file, human genome reference file, regions file
* output : a normal VCF file
* scripts : scripts_cobalt/template_strelka2.sh, scripts_cobalt/launch_strelka2.sh

Then th VCF output file has been filtered on _PASS_ value : keep_pass.sh.

#### Calling quality control

To generate four different charts :

<table>
<thead>
<tr>
  <th>Germline AF distribution</th>
  <th>Somatic Venn</th>
</tr>
</thead>
<tbody>
<tr>
  <td><img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/Calling_quality_control_germline_AF.png" /></td>
  <td><img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/Calling_quality_control_Venn.png" /></td>
</tr>
</tbody>
<thead>
<tr>
  <th>Somatic AF distributions and overlap</th>
  <th>Somatic / Germline overlap</th>
</tr>
</thead>
<tbody>
<tr>
  <td><img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/Calling_quality_control_tumors_overlap.png" /></td>
  <td><img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/Calling_quality_control_tumors_normal_overlap.png" /></td>
</tr>
</tbody>
</table>

* tool : R
* inputs : VCF files
* outputs : SVG files
* scripts : R_controle_qualite/*

#### VCF normalization and annotation

The VCF can be normalized to have a format compatible with annovar.  
* Platypus : mandatory
* Strelka2 : already compatible

###### Normalization

* tool : vt
* inputs : a tumor VCF file
* output : a tumor VCF file
* script : normalize_vcf.sh

###### Annotation

* tool : Annovar
* inputs : a VCF file
* output : a VCF file with extra columns
* script : scripts_cobalt/run_annovar.sh

(In this pipeline, Annovar has been run on somatic VCF only)

#### Tumor coverage

For each patient, somatic calling of tumor1 & 2 give two variant lists with their respective positions and coverage. Germline calling also gives a variant list with its own positions and coverage.        

We need the coverage of these positions in the others samples. For example, we need the coverage in tumor2 at the positions of tumor1 somatic variants. Inversely, we need the coverage in tumor1 at the positions of tumor2 somatic variants.              

We also need tumor1 & 2 coverage at the positions of the germline variants.                    


<img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/coverage.png" />

###### Tumor coverage at the other tumor positions

* tool : Platypus
* inputs : a tumor BAM, a tumor VCF, human genome reference file, regions file
* output : a tumor VCF file
* script : scripts_cobalt/calling_somatic_genotype/platypus_reads_*.sh

###### Tumor coverage at the germline positions

* tool : Platypus
* inputs : a tumor BAM, a normal VCF, human genome reference file, regions file
* output : a tumor VCF file
* script : scripts_cobalt/calling_germline_genotype/platypus_genotype_*.sh

#### Somatic allele-specific copy numbers profiling

###### To get the copy numbers

The script is parallelized by chromosome.  
To split germline VCF by chromosome, use this script : split_vcf_chromosome.sh

* tool : Falcon (R package)
* inputs :
  * args[1] = germline VCF file
  * args[2] = tumor1 sample ID
  * args[3] = tumor2 sample ID
  * args[4] = chromosome
  * args[5] = path to output directory
* outputs :
  * a PDF with segmentation results
  * a PDF with QC
  * a TXT with copy numbers and their standard error
  * a RDA with germline data loaded
  * a RDA with tumor1 copy number computed
  * a RDA with tumor2 copy number computed
* script : marathon/falcon/falcon.R

###### To get the copy numbers, in the other tumor regions with variations

* tool : Falcon (R package)
* inputs :
  * args[1] = tumor1 RDA file
  * args[2] = path to Falcon TXT file containing coordinates (give a path like "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_##/chr6/falcon.patient_##.tumor_placeholder.chr_6.output.txt" so it can take tumor1 and tumor2 TXT file thanks to the placeholder)
  * args[3] = path to output directory
* outputs :
  * a TXT with the computed standard errors
* script : marathon/falcon/falcon_epsilon.R

###### Notes

It is important to note that these Falcon scripts use some custom libraries stored in marathon/libs/.  
Sometimes, these libraries are simple overrides of Falcon with little modifications.

#### Heterogeneity characterization and tree generation

###### SNA pre-clustering and Monte Carlo Markov chain sampling

This step computes all input matrices required by Canopy, and performs SNA pre-clustering, and then a MCMC sampling to give subclones with composition and history.  

This step has been parallelized by number of subclones.

* tool : Canopy (R package)
* inputs :
  * args[1] = path to Falcon patient output directory
  * args[2] = patient ID
  * args[3] = tumor1 ID
  * args[4] = tumor2 ID
  * args[5] = path to tumor1 VCF file
  * args[6] = path to tumor2 VCF file
  * args[7] = path to tumor1 coverage VCF file at tumor2 positions
  * args[8] = path to tumor2 coverage VCF file at tumor1 positions
  * args[9] = path to output directory
  * args[10] = number of subclones to generate
* outputs :
  * a BIC file with bayesian Information Criterion score for this number of subclones
  * a SVG file with optimal number of clusters (non-informative without a comparison with others numbers of clusters)
  * a PDF file with optimal number of subclones (non-informative without a comparison with others numbers of subclones)
  * a PDF file with the likelihood of this number of subclones
  * a RDA file containing the pre-clustering and the MCMC sampling computed for this number of subclones
* script : marathon/canopy/canopy.R

###### Tree generation

* tool : Canopy (R package)
* inputs :
  * args[1] = patient ID
  * args[2] = path to canopy patient output
* outputs :
  * a SVG file with a plot of BIC of each number of subclones (to choose the better BIC)
  * a PDF file with the visual generated tree
  * a TXT file with the composition of each subclones (SNAs and CNAs)
* script : marathon/canopy/canopy_tree.R

###### Events filtering

patient_id            = args[1]
data_path             = args[2]
file_name             = args[3]
input_somatic_VCF_t1  = args[4]
input_somatic_VCF_t2  = args[5]
only_exonic           = args[6]

* tool : Canopy (R package)
* inputs :
  * args[1] = patient ID
  * args[2] = path to canopy patient output
  * args[3] = name of the TXT subclones composition file
  * args[4] = path to the tumor1 VCF file
  * args[5] = path to the tumor2 VCF file
  * args[6] = _1_/_0_. If 1, keep exonic SNAs only
* outputs :
  * a TSV file with filtered events
* script : marathon/canopy/canopy_subclones_composition.R

###### Notes

It is important to note that these Canopy scripts use some custom libraries stored in marathon/libs/.  
Sometimes, these libraries are simple overrides of Canopy with little modifications.

## References

[1] **Assessing intratumor heterogeneity and tracking longitudinal and spatial clonal evolutionary history by next-generation sequencing**  
Jiang, Y., Qiu, Y., Minn, A.J. and Zhang, N.R., 2016.   
Proceedings of the National Academy of Sciences.  
http://www.pnas.org/content/pnas/113/37/E5528.full.pdf      

[2] **Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny**  
Eugene Urrutia Hao Chen Zilu Zhou Nancy R Zhang Yuchao Jiang  
Bioinformatics, bty057 (2018)
https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty057/4838234?redirectedFrom=fulltext
