# Pipeline ITH

## Description

A pipeline to study intratumor heterogeneity (ITH) with Canopy<sup>[1]</sup>.

## General overview

<img src="https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/pipeline_overview.png" style="width:80%; display:block; margin:auto;" />

## Notes

In this documentation, patient ID has been replaced with _##_, and tumors ID has been replaced with _T1_ and _T2_.

## Steps

### Post-alignment

* tool : BWA
* input : a BAM file
* output : a BAM file
* scripts : scripts_cobalt/template_postalt.sh, scripts_cobalt/launch_postalt.sh

### Germline calling

* tool : Platypus
* inputs : a normal BAM file, human genome reference file, regions file
* output : a normal VCF file
* scripts : scripts_cobalt/template_platypus.sh, scripts_cobalt/launch_platypus.sh

Then th VCF output file has been filtered on _PASS_ value : keep_pass.sh.

### Somatic calling

* tool : Strelka2
* inputs : a tumor BAM file, its associated normal BAM file, human genome reference file, regions file
* output : a normal VCF file
* scripts : scripts_cobalt/template_strelka2.sh, scripts_cobalt/launch_strelka2.sh

Then th VCF output file has been filtered on _PASS_ value : keep_pass.sh.

### Calling quality control

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

### VCF normalization and annotation

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

### Tumor coverage at the other tumor positions

* tool : Platypus
* inputs : a tumor BAM, a tumor VCF, human genome reference file, regions file
* output : a tumor VCF file
* script : scripts_cobalt/calling_somatic_genotype/platypus_reads_*.sh

### Tumor coverage at the germline positions

* tool : Platypus
* inputs : a tumor BAM, a normal VCF, human genome reference file, regions file
* output : a tumor VCF file
* script : scripts_cobalt/calling_germline_genotype/platypus_genotype_*.sh

### Somatic allele-specific copy numbers profiling

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

## References

[1] **Assessing intratumor heterogeneity and tracking longitudinal and spatial clonal evolutionary history by next-generation sequencing**  
Jiang, Y., Qiu, Y., Minn, A.J. and Zhang, N.R., 2016.   
Proceedings of the National Academy of Sciences.  
http://www.pnas.org/content/pnas/113/37/E5528.full.pdf
