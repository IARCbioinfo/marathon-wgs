# Pipeline ITH

## Description

A pipeline to study intratumor heterogeneity (ITH) with Canopy<sup>[1]</sup>.

## General overview

![alt text](https://raw.githubusercontent.com/IARCbioinfo/marathon-wgs/master/images/pipeline_overview.png "Pipeline overview")

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

### Somatic calling

* tool : Strelka2
* inputs : a tumor BAM file, its associated normal BAM file, human genome reference file, regions file
* output : a normal VCF file
* scripts : scripts_cobalt/template_strelka2.sh, scripts_cobalt/launch_strelka2.sh

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

## References

[1] **Assessing intratumor heterogeneity and tracking longitudinal and spatial clonal evolutionary history by next-generation sequencing**  
Jiang, Y., Qiu, Y., Minn, A.J. and Zhang, N.R., 2016.   
Proceedings of the National Academy of Sciences.  
http://www.pnas.org/content/pnas/113/37/E5528.full.pdf
