# DisBalance
Data and codes for DisBalance web-server


## Balance transformation
```py
import subprocess


otu_file='./input/D003093/otu_table.xls' 
sample_file='./input/D003093/sample_data.xls.diseaseId.csv'
dlimits=0.0001

Rcmd='Rscript ./bin/run_distal_DBA.R {} {} {}'.format( 
                        otu_file, sample_file, dlimits)
subprocess.run(Rcmd, shell=True)
```
## Model building risk Prediction

refer to [mMAL](http://lab.malab.cn/soft/mAML/download.html)

## The performance of the model
The comparison between the different over-sampling algorithms.  

script: `robust_measure_performance_nested.py`
`output/Compare_overSampling_algorithms/`


Independent datasets test

output dir:  `output/QuinnNestedCV/`


## Biomarker Discovery

```py
import subprocess

lrcoefile='./input/D003093/otu_table.xls.balance.csv.LogisticRegression.hypertuned.coef.newids.csv' 
sbpf='./input/D003093/otu_table.xls.balance.sbp.csv'
MicroPhenoDBfile='./input/MicroPhenoDB.GI.ncbitaxid.xls'
Mesh2diseasef='./input/GMrepoMeshIDdiseaseUniq.xls'
topn=10
outdir='./output'

Rcmdtop='Rscript /opt/services/djangoapp/src/mAML/balance/dashApps/utils/run_topbalance2evidence.R {} {} {} {} {} {}'.format(
    lrcoefile, sbpf, MicroPhenoDBfile,Mesh2diseasef,str(topn),outdir) 

subprocess.run(Rcmdtop, shell=True)
```
