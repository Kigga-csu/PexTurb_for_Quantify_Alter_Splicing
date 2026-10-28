import numpy as np
from multiprocessing import Pool
import multiprocessing
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tpm", "-t", required=True, dest="tpm",
                    help="transcript expression file")
parser.add_argument("--index", "-i", required=True, dest="index",
                    help="Index file for calc_psi input.py")
parser.add_argument("--psi", "-p", required=True, dest="psi",
                    help="psi file for utest.py input")

def division(markers:list, pair_idx, den_idx, fpkm_mat):
    if pair_idx is None:
        pair_cnts = np.zeros(fpkm_mat.shape[1])
    else:
        pair_cnts = np.sum(fpkm_mat[pair_idx], axis=0)
    den_cnts = np.sum(fpkm_mat[den_idx], axis=0)
    psi = pair_cnts / den_cnts
    #print("4",flush=True)
    psi = [str(i) for i in psi]
    return "\t".join(markers+psi)+"\n"

def pair_psi(fpkm_mat, gtf_index_path, sample_id, psi_file):
    with open(gtf_index_path, "r") as f:
        head = f.readline().rstrip().split("\t")[:3] + sample_id
        f_w = open(psi_file, "w")
        f_w.write("\t".join(head)+"\n")
        for line in f:
            line = line.rstrip().split("\t")
            if line[-2] == "":
                pair_idx = None
            else:
                line[-2] = line[-2].split(",")
                print(line[-2])
                pair_idx = np.array(line[-2], dtype="int")
            line[-1] = line[-1].split(",")
            den_idx = np.array(line[-1], dtype="int")
            result = division(line[:3], pair_idx, den_idx, fpkm_mat)
            f_w.write(result)
        f_w.close()
        f.close()
    return

def convert_nan_to_na(file_path):
    # 打开原文件和一个临时文件
    with open(file_path, 'r', encoding='utf-8') as infile, open(file_path + '.tmp', 'w', encoding='utf-8') as outfile:
        for line in infile:
            # 将行分割为列
            columns = line.strip().split('\t')
            # 从第4列开始处理
            for i in range(3, len(columns)):
                if columns[i].lower() == 'nan':
                    columns[i] = 'NA'
            # 将处理后的列重新组合为一行
            outfile.write('\t'.join(columns) + '\n')

    # 用临时文件替换原文件
    import os
    os.replace(file_path + '.tmp', file_path)

if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.tpm, "r") as f:
        sample_id = f.readline().rstrip().split("\t")[1:]
        transcript_mat = []
        for line in f:
            line = line.rstrip().split("\t")
            transcript_mat.append(line[1:])
        transcript_mat = np.array(transcript_mat, dtype="float")
    pair_psi(transcript_mat, args.index, sample_id, args.psi)
    convert_nan_to_na(args.psi)
