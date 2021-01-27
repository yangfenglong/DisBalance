args=commandArgs(T)
lrcoefile = ifelse(args[1]=='',NA, args[1])
sbpf = ifelse(args[2]=='',NA, args[2])

MicroPhenoDBfile = ifelse(args[3]=='',NA, args[3])
Mesh2diseasef= ifelse(args[4]=='',NA, args[4])
topn = ifelse(args[5]=='',20, args[5]) 
outdir=ifelse(args[6]=='','./', args[6]) 
#   otuf='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls'
#   groupmff='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/sample_data.xls.diseaseId.csv'
# lrcoefile='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls.balance.csv.LogisticRegression.hypertuned.coef.newids.csv'
# sbpf='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls.balance.sbp.csv'
# MicroPhenoDBfile='/opt/services/djangoapp/src/mAML/balance/static/MicroPhenoDB/MicroPhenoDB.GI.ncbitaxid.xls'
# Mesh2diseasef='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/GMrepoMeshIDdiseaseUniq.xls'

source("../lib/top20_balance2evidenceHeatmap.R")

if(any(is.na(c(lrcoefile,sbpf)))){
    print('requied files are missing ...')
    q(save='no')
}
Top20balanceMicroPhenoDB_Heatmap(
    lrcoefile=lrcoefile, #'balance.csv.LogisticRegression.hypertuned.coef.csv'
    sbpf=sbpf,
    MicroPhenoDBfile=MicroPhenoDBfile, #"MicroPhenoDB.ncbitaxid.xls", 
    Mesh2diseasef=Mesh2diseasef,
    topn=as.integer(topn),
    outdir=outdir
)
