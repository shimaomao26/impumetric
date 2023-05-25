
**Introduction**

This repo offers `impumetric`, a tool to examine the imputation result. Especially,

* Rsq: The standard Minimac metric.

* EmpRsq: The squared correlation between WGS and imputed datasets. Also known as the dosage $r^2$.

* MARE: MAF-adjusted-residual-error, to quantify the deviation between Rsq and EmpRsq, which indicates some systematic bias in the imputation pipeline.

* $\beta_{imp}$: The regression slope between imputed dosage and true genotype, which indicates the imputed dosages are shrunk to the alternative allele frequency or not.


It also provides some functions to generate plots and calculate other metrics regarding the imputation quality.


Citation:

[*#1* TBD](TBD)



**Usage**


1. Using it as a script


**Vcf file from Minimac4 `--meta` option could be directly parsed.**

* Directly parsing the loovcf

```sh
python impumetric.py --loovcf ../data/chr19.empiricalDose.vcf.gz

head chr19.empiricalDose.vcf.gz.scores.txt
```

The loovcf file contains the true allele and imputed haploid dosage. Rsq, EmpRsq, MARE, and Beta_imp will be calculated on the haploid data.

```
        Rsq     EmpR    MARE    Beta_imp        EmpRsq
chr19:247265:G:T        0.862   0.878   0.197   0.712   0.772
chr19:267213:C:T        0.461   0.739   0.209   0.636   0.546
chr19:277717:G:A        0.828   0.89    0.172   0.877   0.792
```


**In other case, it requires [`bcftools`](https://samtools.github.io/bcftools/) to extract the genotypes/dosages first.** Some examples:

* Extracting the variants*sample table of `GT`

```sh
inpath=$1
outpath=$2
variants=overlapping.155297.idv
sample_file=overlapping.n993.iid

# print
bcftools view \
  --include ID==@$variants \
  --force-samples \
  -S $sample_file \
  -Ou $inpath | bcftools query \
  -f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > $outpath
```

* Extracting all `HDS`

```sh
bcftools query -f'%ID[\t%HDS]\n' $inpath > $outpath
```



*  Providing both WGS and imputed dataset (should have the same variants and samples)

```sh
python impumetric.py --imputed ../data/chr19.HDS.sample.txt --wgs ../data/chr19.GT.sample.txt

head chr19.HDS.sample.txt.scores.txt   
```

WGS dataset is usually not phased, thus, the values will be calculated on the diploid data. (Rsq will be on haploid data when `HDS` is available.)

```   
0       Rsq     EmpR    MARE    Beta_imp        EmpRsq
chr19:247265:G:T        0.862   0.872   0.197   0.709   0.76
chr19:267213:C:T        0.461   0.738   0.212   0.66    0.544
chr19:277717:G:A        0.828   0.89    0.185   0.87    0.792
```

* Providing the imputed dataset only

```sh
python impumetric.py --imputed ../data/chr19.HDS.sample.txt

head chr19.HDS.sample.txt.scores.Rsq.txt
```

If only a imputed dataset is provided, only Rsq could be calculated.

```
0       Rsq
chr19:247265:G:T        0.862
chr19:267213:C:T        0.461
chr19:277717:G:A        0.828
```



2. Using it via `import`

Plotting functions are not parsered. These functions need import.

```py
from scores import *
from plot import *

fig = plot_Beta(newdf, "beta", "EmpRsq", filename=None)
```

**The `ipynb` and `html` files under `src` give some examples**, e.g., matching the variant ID between two datasets, calculating metric values, and plotting.



3. Using it as a python package

It could also be imported as a package and integrated into a pipeline.

```py
import impumetric as imp

imp.plot_scatter()
```

Or:

```py
#dir(imp)

from impumetric.align import *

merge_multi_allelic()
```

PS: import package from any place

```py
import sys
sys.path.append(r"path to impumetric package/src")
import impumetric as imp

imp.plot_scatter("EAS.chr19.info.gz", xcol="EmpRsq", ycol="LooRsq", xname=None, yname=None, 
                 filename="EAS.scatter", contour=True)
```