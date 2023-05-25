# Date 20230314, Copyright at Shi

import argparse
#import os
#import gzip

#from align import *
from scores import *
#from plot import *

parser = argparse.ArgumentParser()

# NOTE: To simplify the workflow, now all vcf processing should be done in bcftools, 
# except directly parsing the loodosage. Please check the scripts in the repo for pre-processing.

#######
# score
# read the processed table or the LooDosage vcf
parser.add_argument("--imputed", default=None, help="the imputed dataset (processed table by bcftools)")
parser.add_argument("--wgs", default=None, help="the wgs dataset (processed table by bcftools)")
parser.add_argument("--loovcf", default=None, help="now only parsing the vcf of loodosage is supported. \
    In this case, --imputed, --wgs are not required.")

# calculate the scores
# NOTE: To simplify the usage, now it will automatically calculate the scores given the specific data input
# parser.add_argument("--scores-imputed", default=None, help="if the haploid data is provided, minimac Rsq will be calculated; \
#     if the diploid data is provided, MaCH Rsq will be calculated")
# parser.add_argument("--scores-wgs", default=None, help="because WGS data is diploid, the correlation between diploid data \
#     will be calculated. Both Pearson r and r^2 will be returned")

#########
# Add-ons
# parser.add_argument("--combine-rsq", default=None, help="Experimental function!! Recalculate scores from batch imputation. \
#     Require a list of score files, with N, AAF, Rsq in each file. This function will return a combined score file. For example, \
#     you could input the Minimac4 *.info files to recalculate AAF and Rsq. Due to the rounding error, the recalculated \
#     scores will be different from submitting all samples to the Imputation Server, especially for ultra rare variants.")
# parser.add_argument("--combine-N", default=None, help="By default, sample size is not in the Minimac4 info file, input it here \
#     in the same order with the --combine-rsq option")

# NOTE: latent functions were not parsed here. To use the latent function, using "import" or directly copy the code.


#####
# run
#####
args = parser.parse_args()

if args.loovcf:
    print("Now it is using the looDosage file")
    true_gt0=pd.read_table(args.loovcf,  comment='#', header=None).iloc[:,9:]
    #true_gt0 = loo_ds
    gt = pd.concat([true_gt0.applymap(lambda x: x.split(":")[0].split("|")[0]), 
                        true_gt0.applymap(lambda x: x.split(":")[0].split("|")[1])], axis=1)
    #lds = pd.read_table(fname1+i+fname2, comment='#', header=None).iloc[:,9:]
    lds = pd.concat([true_gt0.applymap(lambda x: x.split(":")[1].split("|")[0]), 
                    true_gt0.applymap(lambda x: x.split(":")[1].split("|")[1])], axis=1)
    gt = gt.astype("int8") # No missing data should be in the looDosage
    lds = lds.astype("float16")

    Rsq = rsq_minimac(hds=lds)
    EmpR = cor_emp(gt, lds)
    vid, Beta_imp, MARE = mare_beta_hap(args.loovcf)

    df = pd.DataFrame()
    df.index=vid
    df["Rsq"] = Rsq
    df["EmpR"] = EmpR
    df["MARE"] = MARE
    df["Beta_imp"] = Beta_imp
    # if the corr < 0, corr^2 = 0
    df["EmpRsq"] = np.where(df["EmpR"] >=0, df["EmpR"]**2, 0)

    df = df.round(3)

    df.to_csv(args.loovcf+".scores.txt",sep="\t")

if args.imputed:
    if args.loovcf:
        print("Because looDosage file is used, other args will be ignored")

    elif args.wgs:
        print("Trying to calculate EmpRsq between WGS and imputed data")
        gt = pd.read_table(args.wgs,  header=None, index_col=0) # wgs genotype, missing data should be encoded as .
        ds = pd.read_table(args.imputed,  header=None, index_col=0) # imputed ds or hds
        assert len(gt) == len(ds), "wgs and imputed dataset should have the same variants (in the same order)"
        assert gt.shape[1] == ds.shape[1], "wgs and imputed dataset should have the same individuals (in the same order)"

        # if hds is in 1,1 format, calculate Rsq from hds, otherwise, calculate MaCH Rsq from ds
        if type(ds.iloc[0,0]) == str:
            if "," in ds.iloc[0,0]:
                ds_ = pd.concat([ds.applymap(lambda x: x.split(",")[0]), 
                    ds.applymap(lambda x: x.split(",")[1])], axis=1).astype("float16")
                Rsq = rsq_minimac(hds=ds_)
                ds = ds.applymap(lambda x: float(x.split(",")[0]) + float(x.split(",")[1]))
            elif "|" in ds.iloc[0,0]:
                ds_ = pd.concat([ds.applymap(lambda x: x.split("|")[0]), 
                    ds.applymap(lambda x: x.split("|")[1])], axis=1).astype("float16")
                Rsq = rsq_minimac(hds=ds_)
                ds = ds.applymap(lambda x: float(x.split("|")[0]) + float(x.split("|")[1]))
            else:
                raise ValueError("Imputed data format should be any of the follow: 0.224,0.225 0.224/0.225 0.449")
        elif type(ds.iloc[0,0]) == float:
            Rsq = rsq_MaCH(ds)
        else:
            raise ValueError("Format is neither str nor float")

        # gt in 1/1 or 1|1 format
        if type(gt.iloc[0,0]) == str:
            if "/" in gt.iloc[0,0]:
                gt1 = gt.applymap(lambda x: x.split("/")[0]).replace(".", np.nan).astype("float16") #.astype("Int8") # allow missing
                gt2 = gt.applymap(lambda x: x.split("/")[1]).replace(".", np.nan).astype("float16") #.astype("Int8") # Int8 doesn't support np.corrcoef; use float16
                gt = gt1 + gt2
            elif "|" in gt.iloc[0,0]:
                gt1 = gt.applymap(lambda x: x.split("|")[0]).replace(".", np.nan).astype("float16") #.astype("Int8") # allow missing
                gt2 = gt.applymap(lambda x: x.split("|")[1]).replace(".", np.nan).astype("float16") #.astype("Int8")
                gt = gt1 + gt2
            else:
                raise ValueError("WGS data format should be any of the follow: 1|0 1/0 1")
        elif type(gt.iloc[0,0]) == float or type(gt.iloc[0,0]) == int:
            gt = gt.replace(".", np.nan).astype("float16") #.astype("Int8")
        else:
            raise ValueError("Format is neither str nor float")

        # EmpR is always on gt and ds
        #gt = gt.astype("Int8")
        EmpR = cor_emp(gt, ds)
        vid, Beta_imp, MARE = mare_beta_dip(ds, gt)
    
        df = pd.DataFrame()
        df.index=vid
        df["Rsq"] = Rsq
        df["EmpR"] = EmpR
        df["MARE"] = MARE
        df["Beta_imp"] = Beta_imp
        # if the corr < 0, corr^2 = 0
        df["EmpRsq"] = np.where(df["EmpR"] >=0, df["EmpR"]**2, 0)

        df = df.round(3)

        df.to_csv(args.imputed+".scores.txt",sep="\t")

    else:
        print("Trying to calculate Rsq using the imputed data only")
        ds = pd.read_table(args.imputed,  header=None, index_col=0)

        # if hds is in 1,1 format, calculate Rsq from hds, otherwise, calculate MaCH Rsq from ds
        if type(ds.iloc[0,0]) == str:
            if "," in ds.iloc[0,0]:
                ds_ = pd.concat([ds.applymap(lambda x: x.split(",")[0]), 
                    ds.applymap(lambda x: x.split(",")[1])], axis=1).astype("float16")
                Rsq = rsq_minimac(hds=ds_)
                ds = ds.applymap(lambda x: float(x.split(",")[0]) + float(x.split(",")[1]))
            elif "|" in ds.iloc[0,0]:
                ds_ = pd.concat([ds.applymap(lambda x: x.split("|")[0]), 
                    ds.applymap(lambda x: x.split("|")[1])], axis=1).astype("float16")
                Rsq = rsq_minimac(hds=ds_)
                ds = ds.applymap(lambda x: float(x.split("|")[0]) + float(x.split("|")[1]))
            else:
                raise ValueError("Format should be any of the follow: 0.224,0.225 0.224/0.225 0.449")
        elif type(ds.iloc[0,0]) == float:
            Rsq = rsq_MaCH(ds)
        else:
            raise ValueError("Format is neither str nor float")

        df = pd.DataFrame()
        df.index = ds.index
        df["Rsq"] = Rsq

        df = df.round(3)

        df.to_csv(args.imputed+".scores.Rsq.txt",sep="\t")

# if args.combine_rsq:
#     print("Now it is trying to recalculate Rsq from different files containing MAF, Rsq, and N")

print("End of run")