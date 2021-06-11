if __name__ == '__main__': 
    aln_rootDir = '/data/srh/rawdata/'
    SRH_rootDir = '/data/srh/processed_data/SRH_tables/'
    IQtree_rootDir = '/data/srh/processed_data/IQtree/'
    for DirName, subdirList, fileList in os.walk(aln_rootDir):
        for fname in fileList:
            if(fname=="alignment.nex"):
                print(DirName)
                head_DirName, datas = os.path.split(DirName)
                aln_path = os.path.join(DirName,fname)

                # if this file exists, we already did this dataset successfully
                if not os.path.exists(os.path.join(IQtree_rootDir,datas,"MPTMS","Bad")):
                    dat = Nexus.Nexus()
                    dat.read(aln_path) 
                    aln = AlignIO.read(open(aln_path), "nexus")
                    p = Test_aln(aln,datas,dat)
                    TempDir = os.path.join(SRH_rootDir,datas,'Data')
                    if not os.path.exists(TempDir):
                        os.makedirs(TempDir)
                        os.chdir(TempDir)
                    df =pd.DataFrame(p[1:], columns=p[0])
                    T=table(p)
                    T.to_csv(os.path.join(TempDir,'MaxDiv.csv'))
                    df.to_csv(os.path.join(TempDir,'data2.csv'))
                    merged = pd.read_csv(os.path.join(TempDir,'merged.csv'))
                    merged = pd.merge(merged, T, how='outer', on=['dataset','Charset','Test'])
                    merged.drop(['isbad', 'Unnamed: 0'], axis=1, inplace=True)
                    merged.to_csv(os.path.join(TempDir,'merged2.csv'))
                    all_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","All")
                    if not os.path.exists(all_MPTS_path):
                        os.makedirs(all_MPTS_path)
                    good_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","Not_Bad")
                    if not os.path.exists(good_MPTS_path):
                        os.makedirs(good_MPTS_path)
                    bad_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","Bad")
                    if not os.path.exists(bad_MPTS_path):
                        os.makedirs(bad_MPTS_path)
                    all_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","All")
                    if not os.path.exists(all_MPTIS_path):
                        os.makedirs(all_MPTIS_path)
                    good_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","Not_Bad")
                    if not os.path.exists(good_MPTIS_path):
                        os.makedirs(good_MPTIS_path)
                    bad_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","Bad")
                    if not os.path.exists(bad_MPTIS_path):
                        os.makedirs(bad_MPTIS_path)
                    all_MPTMS_path = os.path.join(IQtree_rootDir,datas,'MPTMS','All')
                    if not os.path.exists(all_MPTMS_path):
                        os.makedirs(all_MPTMS_path)
                    good_MPTMS_path = os.path.join(IQtree_rootDir,datas,"MPTMS","Not_Bad")
                    if not os.path.exists(good_MPTMS_path):
                        os.makedirs(good_MPTMS_path)
                    bad_MPTMS_path = os.path.join(IQtree_rootDir,datas,"MPTMS","Bad")
                    if not os.path.exists(bad_MPTMS_path):
                        os.makedirs(bad_MPTMS_path)
                    MPTS_all_file = os.path.join(all_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTS_path,'alignment.nex'))
                    MPTIS_all_file = os.path.join(all_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTIS_path,'alignment.nex'))
                    MPTMS_all_file = os.path.join(all_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTMS_path,'alignment.nex'))
                    MPTS_good_file = os.path.join(good_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTS_path,'alignment.nex'))
                    MPTIS_good_file = os.path.join(good_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTIS_path,'alignment.nex'))
                    MPTMS_good_file = os.path.join(good_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTMS_path,'alignment.nex'))
                    MPTS_bad_file = os.path.join(bad_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTS_path,'alignment.nex'))
                    MPTIS_bad_file = os.path.join(bad_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTIS_path,'alignment.nex'))
                    MPTMS_bad_file = os.path.join(bad_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTMS_path,'alignment.nex'))
                    partition_files(T,aln_path)
