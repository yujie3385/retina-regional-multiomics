# Human retinal regional multiomics

Analysis code for central and peripheral human retina across fetal development and adulthood.
The project compares regional transcriptional programs, functional modules, candidate chromatin regulatory associations and predicted transcription factor perturbation responses.

## Repository structure

```text
scripts/snRNA/   RNA processing, pseudobulk DESeq2, trajectories, GSEA, AUCell, validation and virtual knockout
scripts/snATAC/  ATAC processing, chromVAR, candidate TF–peak–gene networks, footprints and locus plots
docs/           Environment and reproducibility notes
script_manifest.csv  Published script paths and SHA256 checksums
```

Scripts retain the numbered order of the final analysis collections. R notebooks use an R kernel. The two Python scripts assemble marker-peak figures from previously generated panels.

## Data

Primary RNA resources are available through the [Human Retina Cell Atlas](https://rchenlab.github.io/resources/human-atlas.html).
ATAC fragments and sample information are available under [GSE265801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265801).
Independent validation uses GSE142526, GSE277326 and GSE155288.

RNA count matrices, cell annotations and ATAC peak matrices will be linked here when their archive location is available. They are not included in this repository. The publication release is intended to support a matrix-based starting point; the supplied scripts also retain upstream processing code for provenance.

Matrix columns must match annotation cell identifiers. Annotation tables must identify cell type, donor, sample, age, developmental stage and retinal region. ATAC inputs also require hg38 peak coordinates. Coverage tracks and motif footprints require fragment files and their tabix indexes; a peak matrix alone cannot reproduce these panels.

## Paths and execution

Personal usernames, original server directory prefixes and local drive paths have been removed. Notebook outputs and execution metadata have been cleared before publication.

The following neutral prefixes are placeholders, not existing directories:

| Prefix | Intended contents |
| --- | --- |
| `/project/work/` | Analysis working directories and intermediate results |
| `/project/resources/` | Reference objects, metadata and external resources |
| `/project/local/` | Inputs previously located on a local workstation |
| `/project/home/` | Auxiliary resources previously under a home directory |

Configure the input and output paths near the beginning of each script for your environment. Some scripts contain additional path literals. Relative paths and original intermediate file basenames are retained where possible. Local workstation paths have been reduced to neutral filenames and must be configured explicitly. This release does not claim that redaction alone makes the collection executable without input preparation.

Heavy computation was performed in a Linux R environment. Run analysis scripts on a suitable compute server, not by loading all single-cell objects on a typical workstation. Scripts write to their configured output directories; use a separate working directory when reproducing results.

## Analysis order

1. Prepare annotated RNA and ATAC inputs and review sample eligibility.
2. Process fetal and adult datasets separately and summarize quality and cell composition.
3. Aggregate RNA counts by cell type, donor, stage and region; perform paired regional DESeq2 comparisons.
4. Classify DM, DAL and AAM patterns, compute AUCell scores and assess independent datasets.
5. Run preranked GSEA and connect leading-edge genes to curated functional modules.
6. Compute regional chromVAR differences and assemble motif-containing peak-to-gene associations and candidate networks.
7. Generate locus tracks and footprints, and fit virtual knockout models in the appropriate RNA groups.
8. Generate figure panels from the corresponding analysis results.

Within snRNA, scripts 01–08 cover processing and overview, 09–16 differential expression and classification, 17–28 signatures and function, 29–36 virtual perturbation, and 37–54 validation and expression displays. Within snATAC, script 00 covers object construction, 01–05 overview and marker panels, 06–15 motif and network analyses, 16–18 footprints, and 19–21 integrated displays and tracks. Numbering describes the collection; some later scripts consume results from earlier modules across modalities.

## Interpretation and final analysis conventions

- Central combines foveal and macular labels; Peripheral retains peripheral tissue.
- RNA stages are Early at 70–91 days, Mid at 100–116 days, Late at 137–165 days, Young at 21–32 years and Old at 67–86 years, subject to cell-type and paired-sample eligibility.
- DM means developmentally maintained; DAL means developmental aging-lost; AAM means adult-associated and maintained. Historical internal class strings in the scripts have been retained to preserve compatibility with saved result tables.
- The final RNA signature displays exclude the unpaired 70-year-old sample. Final adult regional ATAC comparisons and locus tracks use paired donors aged 21 and 67 years.
- Candidate motif and peak-to-gene associations do not establish TF binding or causality. Virtual knockout results describe predicted differential regulation.
- Regional virtual knockout totals use pooled adult networks fitted independently within each region. Young and Old models are separate analyses; regional totals are not a union of stage-specific significant-gene lists.

## Reproducibility status

This initial release preserves the supplied final script collections with path redaction and notebook-output removal. It has not been rerun end to end after redaction.
Some saved script settings predate the final figure choices. In particular, verify the final 70-year exclusion in AUCell scripts, restrict adult locus-track inputs to 21 and 67 years, and use stage grouping for the revised supplementary MG expression panel. These choices are documented here rather than silently changing analysis code without its input objects.

The current collection includes pooled KO result-reading and plotting code; the pooled network-generation step must be supplied or reconstructed from its original execution record before that panel is reproducible from raw RNA counts. The exact sample scope of the original cached shared-peak calling remains unverified. The object-construction script can load an existing peak cache, and its fallback peak-calling section requires the fragment-path list to be defined first. Curated pathway-to-module mappings and other precomputed resources must also be supplied where a script reads them.

## Software

See [docs/software.md](docs/software.md) for recorded versions. No automated installation is performed. Large data, local credentials, generated outputs and private manuscript files are excluded.

## Citation and licensing

Manuscript citation and persistent data/code archive links will be added when available. No software license has been selected for this initial deposit. Original data and third-party resources retain their respective terms.
