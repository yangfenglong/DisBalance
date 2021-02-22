import os, sys
file_dir='/opt/services/djangoapp/src/mAML/balance'

sys.path.append(file_dir) #for dashApps import

utilspath = file_dir + '/dashApps/utils'
sys.path.append(utilspath) 
from dashApps.utils.components import * # otu2balance
from sklearn_pipeline import *


def robust_measure_performance(otu_file, sample_file, LRmodel, LRgrid,
                                prevalence=0.02, dlimits=0.0001,
                                rebalance="SVMSMOTE",transpose=True,to_relative=False,
                                scorer='accuracy', inner_n_splits=10, inner_test_size=0.3, outer_cv=5, inner_cv=5,
                                n_splits=10,test_size=0.33 # StratifiedKFold for independent validation
                                ):
    X, Y =  get_data(otu_file,sample_file)
    if transpose:
        X=X.T
    X.index=Y.index #for Quin test
    if to_relative :
        X = (X.T / X.sum(axis=1)).T 
    Xall = X.loc[X.sum(axis=1)!=0 , X.sum(axis=0)!=0]    
    # filter low prevalence
    if prevalence !=0:
        label_counts = Counter(Y.iloc[:,0])
        within_class_prevalence =[np.sum(Xall[Y.iloc[:,0].values==k]>0, axis=0)/v for k,v in label_counts.items()] 
        Xall = Xall.loc[:, pd.DataFrame(within_class_prevalence).max() > prevalence] #filter low within class prevalence features
        Xall = Xall.loc[Xall.sum(axis=1)!=0 ,:] # filter 0 
    Yall = Y.loc[Xall.index,:]
    if Xall.shape[1] < 6: #特征过滤太多了，不能计算balance
        exit('Error: number of features are less than 5, set a lower prevalence threshold')

    skf = StratifiedKFold(n_splits=n_splits)
    inner_skf =StratifiedKFold(n_splits=inner_n_splits)
    # inner_skf = StratifiedShuffleSplit(n_splits=inner_n_splits, test_size=inner_test_size, random_state=0)
    # skf = StratifiedShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=0)
    # skf = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=0)
    
    X_filename=otu_file[:-4]    
    Xname= X_filename + ".filter_{}.csv".format(prevalence)
    Yname=X_filename + ".filter_{}.mf.csv".format(prevalence)
    Xall.T.to_csv(Xname)
    Yall.to_csv(Yname)
    
    # transform to balances
    Rcmd='Rscript /opt/services/djangoapp/src/mAML/balance/dashApps/utils/run_distal_DBA.R {} {} {}'.format(
            Xname, Yname, dlimits) 
    os.system(Rcmd)
    Xname=Xname[:-4]
    Balance_filename=Xname+'.balance.csv'
    SBP_filename=Xname+'.sbp.csv'
    # sbp=pd.read_csv(SBP_filename, index_col=0, header=0)
    Xall, Yall = get_data(Balance_filename,Yname)
    outdir = os.path.dirname(Balance_filename)
    Xfilename = os.path.basename(Balance_filename)
    target_names=[re.sub(r" +", "_",name) for name in np.unique(Yall.values)]
    skpimb = SklearnPipeline(Xfilename, Xall, Yall, outdir=outdir)
    skpimb.Y = LabelEncoder().fit_transform(skpimb.Y.values.ravel())
    if rebalance != 'None':   
        skpimb.over_sampling(method=rebalance) 
        skpimb.Y = [target_names[i] for i in skpimb.Y]
        Ydf=pd.DataFrame(skpimb.Y)
        Ydf.index=skpimb.X.index
        skpimb.logger.handlers = []
        Xall, Yall = skpimb.X, Ydf
    for nsplit, (train_index, test_index) in zip(range(n_splits),skf.split(Xall, Yall)):
        X_filename=Balance_filename[:-4] + str(nsplit)
        X, X_test = Xall.iloc[train_index,:], Xall.iloc[test_index,:]
        Y, Y_test = Yall.iloc[train_index,:], Yall.iloc[test_index,:]
        Xtrainf= X_filename + ".train.csv" 
        Xtestname= X_filename + ".test.csv"  
        X.to_csv(Xtrainf)
        X_test.to_csv(Xtestname)
        Y = Y.loc[X.index,:] #sort and select samples after feature selection
        Yname=X_filename + ".train.mf.csv"
        Y.to_csv(Yname)  
        Ytestname=X_filename + ".test.mf.csv"        
        Y_test.to_csv(Ytestname)       
        # otu2balance_main
        # Y = LabelEncoder().fit_transform(Y.values.ravel())
        # plot_tsne(X, Y=Y,targets=target_names, filename=otu_file[:-4])
        log_file = Xtrainf[:-4] + ".sklearnPipeline.log"
        
        # train the model
        Xfilename=os.path.basename(Xtrainf)
        skp = SklearnPipeline(Xfilename, X, Y, log=log_file, outdir=outdir) 
        skp.Y = LabelEncoder().fit_transform(skp.Y.values.ravel()) 
        # plot_tsne(skp.X, Y=skp.Y,targets=skp.target_names, filename=skp.filename)
        # if rebalance != 'None':   
        #     skp.over_sampling(method=rebalance)
        skp.hypertune_classifierCV(
            eval(LRmodel), eval(LRgrid), clf_name="LogisticRegression", 
            pltcfm=False, scorer=scorer,inner_cv=inner_skf,n_jobs=1) #outer_cv=outer_cv,
        # test the performance
        y_pred=skp.best_estimator_.predict(X_test)
        cfm_html = X_filename + '.test.LogisticRegression.PyCM_report'
        target_names=[re.sub(r" +", "_",name) for name in np.unique(Y_test.values)]
        Y_test = LabelEncoder().fit_transform(Y_test.values.ravel())
        dic=dict(zip(np.unique(Y_test),target_names))
        actual_vector = [dic[i] for i in Y_test]
        predict_vector = [dic[i] for i in y_pred]
        cm = ConfusionMatrix(actual_vector=actual_vector, predict_vector=predict_vector)   # pycm
        cm.save_html(cfm_html) # cross prediction result
        skp.logger.handlers = []
        #logger.removeHandler(logger.handlers)



LRmodel="LogisticRegression(penalty='elasticnet', l1_ratio=0.15, solver='saga', multi_class='auto',random_state=0)"
LRgrid='dict(C=np.logspace(-4, 4, 3))' # 0.001 to 5, np.logspace(-4, 4, 3))'


# test for one dataset
#otu_file='balance/measure_performance/D003093/otu_table.xls'
#sample_file='balance/measure_performance/D003093/sample_data.xls.diseaseId.csv'
# robust_measure_performance(otu_file, sample_file, LRmodel, LRgrid,
#                                 prevalence=0.2, dlimits=0.0001,
#                                 rebalance="SVMSMOTE",transpose=True,
#                                 scorer='accuracy', outer_cv=10,inner_cv=5,
#                                 n_splits=10)


# test independent datasets from Quinn's distalDBA paper
quinnlist = './quinn.list'
with open(quinnlist) as lst:
    otufs=[otuf.strip() for otuf in lst.readlines()]
    mapfs=[otuf.replace('x.csv','y.csv') for otuf in otufs]

# 3a 3b 3c 3d
for otu_file, sample_file in zip(otufs,mapfs):
    robust_measure_performance(otu_file, sample_file, LRmodel, LRgrid,
                                prevalence=0.02, dlimits=0.0001,
                                rebalance="SVMSMOTE",transpose=False,to_relative=True,
                                scorer='accuracy', #inner_cv=3,#outer_cv=10,roc_auc 
                                n_splits=10, #test_size=0.05,
                                inner_n_splits=5, #inner_test_size=0.3
                                )
# 
# others[otufs[i] for i in [0,3]],[mapfs[i] for i in [0,3]]
for otu_file, sample_file in zip(otufs[3:4],mapfs[3:4]):
    robust_measure_performance(otu_file, sample_file, LRmodel, LRgrid,
                                prevalence=0.02, dlimits=0.0001,
                                rebalance="SVMSMOTE",transpose=False,to_relative=True,
                                scorer='accuracy', #inner_cv=3,#outer_cv=10,roc_auc 
                                n_splits=10, #test_size=0.1,
                                inner_n_splits=10, #inner_test_size=0.3
                                )


###  comparison over-sampling methods #######

D003093_dir='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/DisBalance/DisBalance/output/Compare_overSampling_algorithms/'
Balance_filename=D003093_dir + 'otu_table.xls.filter_0.02_prevalence.balance.csv'
sample_mff=D003093_dir + 'otu_table.xls.filter_0.02_prevalence.mf.csv'
Xfilename=os.path.basename(Balance_filename)

OVERSAMPLES = ['None','ADASYN', 'KMeansSMOTE', 'SMOTE', 'SVMSMOTE', 'BorderlineSMOTE', 'RandomOverSampler']

X, Y =  get_data(Balance_filename, sample_mff)
prevalence=0.0; dlimits=0.0001
scorer='accuracy'; outer_cv=10; inner_cv=5

for i in OVERSAMPLES:
   outdir='D003093_dir'+i
   log_file = outdir + "/sklearnPipeline.log"
   skp = SklearnPipeline(Xfilename, X, Y, log=log_file, outdir=outdir) 
   skp.filter_low_prevalence_features(prevalence=0) #percent of prevalence 
   plot_tsne(skp.X, Y=skp.Y,targets=skp.target_names, filename=skp.filename)
   skp.over_sampling(method=i)
   skp.hypertune_classifier(eval(LRmodel), eval(LRgrid), clf_name="LogisticRegression", 
       pltcfm=True, scorer=scorer, outer_cv=outer_cv, inner_cv=inner_cv,n_jobs=1)

            
