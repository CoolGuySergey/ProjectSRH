def init_partition_files(partition_file):
    with open(partition_file,'w') as f:
        f.writelines("#nexus\n")
        f.writelines("begin sets;\n")
    return

def end_partition_files(partition_file):    
    with open(partition_file,'a') as f:
        f.writelines("end;")
    return


def partition_files(T,aln_path):
    init_partition_files(MPTS_all_file)
    init_partition_files(MPTIS_all_file)
    init_partition_files(MPTMS_all_file)
    init_partition_files(MPTS_good_file)
    init_partition_files(MPTIS_good_file)
    init_partition_files(MPTMS_good_file)
    init_partition_files(MPTS_bad_file)
    init_partition_files(MPTIS_bad_file)
    init_partition_files(MPTMS_bad_file)
    
    df = T.groupby(['Test']).get_group('MPTS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTS_good_file,'a') as good_MPTS:
            good_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTS_all_file,'a') as all_MPTS:
            all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTS_bad_file,'a') as bad_MPTS:
            bad_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTS_all_file,'a') as all_MPTS:
            all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTIS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTIS_good_file,'a') as good_MPTIS:
            good_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTIS_all_file,'a') as all_MPTIS:
            all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTIS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTIS_bad_file,'a') as bad_MPTIS:
            bad_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTIS_all_file,'a') as all_MPTIS:
            all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTMS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTMS_good_file,'a') as good_MPTMS:
            good_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTMS_all_file,'a') as all_MPTMS:
            all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTMS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTMS_bad_file,'a') as bad_MPTMS:
            bad_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTMS_all_file,'a') as all_MPTMS:
            all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    
    end_partition_files(MPTS_all_file)
    end_partition_files(MPTIS_all_file)
    end_partition_files(MPTMS_all_file)
    end_partition_files(MPTS_good_file)
    end_partition_files(MPTIS_good_file)
    end_partition_files(MPTMS_good_file)
    end_partition_files(MPTS_bad_file)
    end_partition_files(MPTIS_bad_file)
    end_partition_files(MPTMS_bad_file)
    
    return
