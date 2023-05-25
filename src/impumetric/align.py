# Date 20230314, Copyright at Shi

# Align the variant ID between the imputed and WGS dataset.

# This script is to handle the major-minor allele swap, inverted region between hg19 and hg38, indel normalization etc, 
# after Liftover, and to annotate the variant ID with CHR:POS:REF:ALT (not rsid).

import numpy as np
import pandas as pd
# import gzip
# import os



def merge_multi_allelic(df1,df2,pos1,pos2,ref1,ref2,alt1,alt2,save,vtype="snv",to_memory=True,verbose=True):
    '''This function is to merge multi-allelic site with exact match between ref and alt, in case there is both 
    major-minor allele swap and +- chain flipping between hg19 and hg38.

    All alleles should be on + chain of that build, and left-normalized.

    Duplicated records should be removed before this function.

    One chromosome should be used.
    
    Args:
        df1, df2: merge df2 to df1
        pos1, pos2: colname for position
        ref1, ref2: colname for the ref allele
        alt1, alt2: colname for the alt allele
        type:
            snv: match ref==ref, alt==alt and ref==alt, alt==ref (a major-minor allele flip between genome builds)
            indel: match ref==ref, alt==alt. Can't determine it is insertion or deletion when observing a ref==alt, alt==ref.
        save: save the matched table to txt like "save_p4p5.txt".
        to_memory: if True, return the matched table to memory. Default to True.
        verbose: print more information

    Return:
        p1: variants with ref==ref, alt==alt (pattern 1)
        p2: variants with ref==alt, alt==ref (pattern 2, usually because of the major-minor allele swap between hg19 and hg38)
        p3: not matched

        p4p5p6 are just for checking and should not be used in production
        p4: variants with ref==ref, alt!=alt (pattern 4, the alt alleles are different between datasets)
        p5: variants with ref!=ref, alt==alt (pattern 5, the ref alleles are different between datasets, usually the ref allele 
        is also the minor allele)
        p6: variants with ref!=ref, alt!=alt (pattern 6, these variants should be removed)
    '''

    if verbose:
        print("{} variants in df1, {} variants in df2".format(len(df1),len(df2)))

    # ref=ref, alt=alt, p1
    _df_1 = df1.reset_index().merge(df2,how="inner",left_on=[pos1,ref1,alt1],right_on=[pos2,ref2,alt2]).set_index("index").copy()
    _df_1["pattern"]=1
    print("{} unique variant with ref==ref, alt==alt (pattern 1)".format(len(np.unique(_df_1.index))))
    if verbose:
        print("{} records matched. Generally, number of records should equal the number of variants".format(len(_df_1)))

    # ref=alt, alt=ref (major minor allele switch), p2
    _df_2 = df1.reset_index().merge(df2,how="inner",left_on=[pos1,ref1,alt1],right_on=[pos2,alt2,ref2]).set_index("index").copy()
    _df_2["pattern"]=2
    print("{} unique variant with ref==alt, alt==ref (pattern 2)".format(len(np.unique(_df_2.index))))
    if verbose:
        print("{} records matched. Generally, number of records should equal the number of variants".format(len(_df_2)))
        print("however, liftover may lift different variants to the same position. Please take caution if variants != records.")

    # p1+p2
    #p1p2_index=_df_1.index.union(_df_2.index) 
    p1p2_index=_df_1.index.symmetric_difference(_df_2.index) 
    if vtype=="snv":
        assert len(p1p2_index)==len(_df_1)+len(_df_2), "some variants show both pattern 1 and pattern 2, only duplicate records could cause this"
    if vtype=="indel":
        print("{} multi-allelic sites show both pattern 1 and pattern 2".format(len(_df_1)+len(_df_2)-len(p1p2_index)))
    # note: index.union will remove intersection between _df_1 and _df_2, 
    # _df_1 and _df_2 do not have the same SNV, but may have indel (with major minor allele swap between p1 and p2)
    # symmetric_difference will remove these indel, aka, multi-allelic indel at the same position will be removed

    # not matched, p3
    _df_3 = df1.loc[df1.index.difference(p1p2_index),:].copy()
    print("{} variants not matched".format(len(_df_3)))
    #assert len(np.unique(p1p2_index))+len(_df_3)==len(df1), "length of p1p2p3 doesn't match df1, something failed" 

        # Deprecated. Matching multi-allelic sites
        # # ref=ref, alt!=alt (different alt allele)
        # # will cause multiple merge problem
        # _df_4=_df_3.reset_index().merge(df2,how="inner",left_on=[pos1,ref1],right_on=[pos2,ref2]).set_index("index").copy()
        # _df_4["pattern"]=4
        # print("{} variants with ref==ref, alt!=alt (pattern 4)".format(len(np.unique(_df_4.index))))
        # print("with multiple merge, {} variants in df4".format(len(_df_4)))

        # # ref!=ref, alt=alt (major minor allele switch + different alt allele)
        # # will cause multiple merge problem
        # _df_5=_df_3.reset_index().merge(df2,how="inner",left_on=[pos1,alt1],right_on=[pos2,alt2]).set_index("index").copy()
        # _df_5["pattern"]=5
        # print("{} variants with ref!=ref, alt==alt (pattern 5)".format(len(np.unique(_df_5.index))))
        # print("with multiple merge, {} variants in df5".format(len(_df_5)))

        # p4p5_index=_df_4.index.union(_df_5.index) # note: index.union will remove intersection between _df_4 and _df_5, 
        # # but not duplicates within _df_4 or _df_5

        # # not matched, p6
        # # ref!=ref, alt!=alt (inverted region (+- chain between hg19 and hg38), or other problem)
        # _df_6 = _df_3.loc[_df_3.index.difference(p4p5_index),:].copy()

        # print("{} variants not matched".format(len(_df_6)))
        # assert len(np.unique(p4p5_index))+len(_df_6)==len(_df_3),"length of p4p5p6 doesn't match p3, something failed"

    if save:
        print("{}_p1p2.txt.gz contains overlapping variants to calculate dosage r^2".format(save))
        # if verbose:
        #     print("! for indel, this file likes to contain AAC/A, A/AAC, which can't be determined it is insert or deletion")
        #     print("! by checking MAF concordance in p1 and p2, it would provide some information !")
        _ = pd.concat([_df_1,_df_2]).loc[p1p2_index, :]
        _.to_csv(save+"_p1p2.txt.gz",sep="\t",index=None)

        #if vtype=="snv":
        if vtype=="snv" or vtype=="indel":
            #print("FOR SNV")
            print("{}_p1p2p3.txt contains overlapping + df1 unique variants to calculate the coverage".format(save))
            #print("! this file do not contain any duplicates")
            pd.concat([_,_df_3]).to_csv(save+"_p1p2p3.txt.gz",sep="\t",index=None,na_rep="NA")
        
        # if vtype=="indel":
        #     #print("FOR INDEL")
        #     _ = df1.loc[df1.index.difference(_df_1.index),:]
        #     print("{}_p1p3.txt contains overlapping + df1 unique variants to calculate the coverage".format(save))
        #     #print("! this file do not contain any duplicates")
        #     pd.concat([_df_1,_]).to_csv(save+"_p1p3.txt",sep="\t",index=None,na_rep="NA")

        # Deprecated. Matching multi-allelic sites
        # print("{}_p4p5.txt contains the multi-allelic site in which alt is different, but ref is same".format(save))
        # print("{}_p4p5.txt did NOT remove multiple merge, j_pos, j_ref, j_alt may have identical values".format(save))
        # print("! for snv, this file is to check the different alt allele imputed in TOPMed and EAS")
        # pd.concat([_df_4,_df_5]).to_csv(save+"_p4p5.txt",sep="\t",index=None)

        # if len(_df_6)>0:
        #     print("{}_p6.txt is df1 only variant. No variant in df2, or in df2 ref!=ref, alt!=alt".format(save))
        #     print("! for indel, this file may lack something... like AAC/A, A/AAC, which means insert in jewel and deletion in topmed")
        #     print("! for indel, in other words, don't use this file for un-imputable indel")
        #     _df_6.to_csv(save+"_p6.txt",sep="\t",index=None)

        #     print("{}_p1p2p4p5p6 should only be used to calculate the coverage".format(save))
        #     #_=pd.concat([_df_1,_df_2,
        #     #_df_4[~_df_4.index.duplicated(keep='first')],
        #     #_df_5[~_df_5.index.duplicated(keep='first')],
        #     #_df_6])
        #     print("! this file contains multiple merge, ambiguous indel")
        #     print("! this file is for calculating the coverage")
        #     print("! before doing that, use AAF to filter the wrongly matched variants")
        #     _=pd.concat([_df_1,_df_2,
        #     _df_4,
        #     _df_5,
        #     _df_6])
        #     _.to_csv(save+"_p1p2p4p5p6.txt",sep="\t",index=None,na_rep="NA")
        # else:
        #     print("{}_p1p2p4p5 should only be used to calculate the coverage".format(save))
        #     #_=pd.concat([_df_1,_df_2,
        #     #_df_4[~_df_4.index.duplicated(keep='first')],
        #     #_df_5[~_df_5.index.duplicated(keep='first')]])
        #     print("! this file contains multiple merge, ambiguous indel")
        #     print("! this file is for calculating the coverage")
        #     print("! before doing that, use AAF to filter the wrongly matched variants")
        #     _=pd.concat([_df_1,_df_2,
        #     _df_4,
        #     _df_5])
        #     _.to_csv(save+"_p1p2p4p5.txt",sep="\t",index=None)

    print("=== final report ===")
    print("{} variants in df1".format(len(df1)))
    #print("{} variants in {}_p1p2p4p5p6.txt".format(len(_),save))
    #print("crude coverage: {} / {}".format(len(np.unique(p1p2_index))+len(_df_2), len(_)))
    print("{} variants in df1 matched to df2 with ref==ref, alt==alt".format(len(_df_1)))
    print("{} variants in df1 matched to df2 with ref==alt, alt==ref".format(len(_df_2)))
    #if vtype=="snv":
    if vtype=="snv" or vtype=="indel":
        print("crude coverage: {} / {}".format(len(np.unique(p1p2_index)), len(df1)))
    # if vtype=="indel":
    #     print("crude coverage: {} / {}".format(len(_df_1), len(df1)))
    #print("crude different alt: {} / {}".format(len(np.unique(p4p5_index))+len(_df_2), len(_)))

    if to_memory:
        return _
    else:
        return True