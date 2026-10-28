import collections
import numpy as np
import scipy.stats as stats
import argparse
from multiprocessing.dummy import Pool
from gtf_index import read_gtf
from statsmodels.stats.multitest import fdrcorrection

parser = argparse.ArgumentParser()
parser.add_argument("--psi", "-p", required=True, dest="psi",
                    help="psi file from the output of psi.py")
parser.add_argument("--output", "-o", required=True, dest="output",
                    help="the path of output file")
parser.add_argument("--anno", "-a", required=True, dest="anno",
                    help="the annotation file of samples")
parser.add_argument("--min_group_cnt", "-m", dest="min_group_cnt",
                    help="min_gruop_cnt", default=1, type=int)


def group(anno_file: str):
    group2samples = collections.defaultdict(list)
    with open(anno_file, "r") as f:
        i = 0
        for line in f:
            i+=1
            #print(i, line)
            line = line.rstrip().split()
            group2samples[line[1]].append(line[0])
        f.close()
    #print(len(group2samples))
    return list(group2samples.values())[0], list(group2samples.values())[1]

def emp_logit(x, eps=1e-3):
    x = np.clip(x, eps, 1 - eps)
    return np.log(x / (1 - x))

def utest(psi_file, group1, group2, output):
    f_result = open(output, "w")
    f_result.write("chr:exon_up+exon_down\tgene_id\tgene_name\tdelta_psi\tpvalue\tfdr\n")
    result = []
    pvalues = []
    with open(psi_file, "r") as f:
        head = f.readline()
        samples = head.rstrip().split("\t")[3:]
        group1_idx = []
        group2_idx = []
        for i in range(len(samples)):
            if samples[i] in group1:
                group1_idx.append(i)
            elif samples[i] in group2:
                group2_idx.append(i)
        group1_idx = np.array(group1_idx).astype("int")
        group2_idx = np.array(group2_idx).astype("int")
        #print(group1_idx)
        #print(group2_idx)
        #print(len(group1)+len(group2))

        for line in f:
            #pool.apply_async(utest_single_line, (line))
            event = line.rstrip().split("\t")[:3]
            psi_strings = line.rstrip().split("\t")[3:]
            psi_strings = [float(x) if x != 'NA' else np.nan for x in psi_strings]
            psi = np.array(psi_strings, dtype="float")
            print(len(psi))
            group1_psi = psi[group1_idx]
            group2_psi = psi[group2_idx]
            #nans = np.isnan(a)
            group1_psi = np.delete(group1_psi, np.where(np.isnan(group1_psi)))
            group2_psi = np.delete(group2_psi, np.where(np.isnan(group2_psi)))
            if len(group1_psi) < args.min_group_cnt or len(group2_psi) < args.min_group_cnt:
                continue
            delta_psi = np.mean(group1_psi) - np.mean(group2_psi)
            # Welch (unequal-variance) t-test on empirical-logit(PSI): continuous p-value
            la, lb = emp_logit(group1_psi), emp_logit(group2_psi)
            if np.ptp(la) == 0 and np.ptp(lb) == 0:
                pvalue = 1.0
            else:
                pvalue = stats.ttest_ind(la, lb, equal_var=False).pvalue
            if not np.isfinite(pvalue):
                pvalue = 1.0
            pvalues.append(pvalue)
            result.append(event + [str(delta_psi), str(pvalue)])
        #print(pvalues)
        _, fdrs = fdrcorrection(pvalues)
        #print(fdrs)
        for i in range(len(result)):
            f_result.write("\t".join(result[i] + [str(fdrs[i])]) + "\n")
    f_result.close()


if __name__ == '__main__':
    args = parser.parse_args()
    group1, group2 = group(args.anno)
    print("group1")
    print(group1)
    print("group2")
    print(group2)
    utest(args.psi, group1, group2, args.output)

