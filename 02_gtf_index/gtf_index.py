import collections
from multiprocessing.dummy import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tpm", "-t", required=True, dest="tpm",
                    help="transcript expression file")
parser.add_argument("--gtf", "-g", required=True, dest="gtf",
                    help="gene annotation file")
parser.add_argument("--index", "-i", required=True, dest="index",
                    help="Index file for calc_psi input.py")


def attributes2dict(attr: list):
    attri_dict = collections.defaultdict()
    for i in attr:
        i = i.split(" ")
        attri_dict[i[0]] = i[1][1:-1]
    return attri_dict

class gtf_dict:
    def __init__(self, gene2transcripts=None,
                 transcript2exons=None,
                 gene2chr=None,
                 geneid2name=None,
                 exon2transcripts=None,
                 site2transcripts=None):
        self.gene2transcripts = gene2transcripts
        self.transcript2exons = transcript2exons
        self.gene2chr = gene2chr
        self.geneid2name = geneid2name
        self.exon2transcripts = exon2transcripts
        self.site2transcripts = site2transcripts

def read_gtf(gtf_path: str):
    gene2transcripts = collections.defaultdict(list)
    transcript2exons = collections.defaultdict(list)
    gene2chr = collections.defaultdict(str)
    geneid2name = collections.defaultdict(str)
    exon2transcripts = collections.defaultdict(list)
    site2transcripts = collections.defaultdict(list)
    chr = range(1, 23)
    chr = [str(i) for i in chr] + ["X", "Y"]
    chr = ["chr"+i for i in chr] + chr
    print(chr)
    with open(gtf_path, "r") as f:
        for line in f:
            line = line.rstrip().split("\t")
            if len(line) != 9 or line[0] not in chr:
                continue
            if line[2] != "exon" and line[2] != "transcript" and line[2] != "gene":
                continue
            line[8] = line[8].split("; ")
            line[8] = attributes2dict(line[8])
            if line[2] == "gene":
                gene = line[8]["gene_id"]
                gene_name = line[8].get("gene_name", gene)
                gene2chr[gene] = line[0]
                geneid2name[gene] = gene_name
            if line[2] == "transcript":
                gene = line[8]["gene_id"]
                transcript = line[8]["transcript_id"]
                gene2transcripts[gene].append(transcript)
            if line[2] == "exon":
                transcript = line[8]["transcript_id"]
                gene = line[8]["gene_id"]
                exon = [line[3], line[4]]
                transcript2exons[transcript].append(exon)
                exon2transcripts[(gene, line[3], line[4])].append(transcript)
                site2transcripts[(gene, "l", line[3])].append(transcript)
                site2transcripts[(gene, "r", line[4])].append(transcript)
    for key, value in transcript2exons.items():
        transcript2exons[key] = sorted(value, key= lambda x: x[0], reverse=False)
    return gtf_dict(gene2transcripts,
                    transcript2exons,
                    gene2chr,
                    geneid2name,
                    exon2transcripts,
                    site2transcripts)

def pair_t(gene2transcripts: dict, transcript2exons: dict):
    def pair2trans(transcripts, transcript2exons):
        pair2trans = collections.defaultdict(list)
        for transcript in transcripts:
            exons = transcript2exons[transcript]
            if len(exons) == 1:
                #print(transcript)
                continue
            for i in range(len(exons)):
                for j in range(i+1, len(exons)):
                    exon1 = "-".join(exons[i])
                    exon2 = "-".join(exons[j])
                    pair2trans[exon1 + "+" + exon2].append(transcript)
        return pair2trans

    pool = Pool(processes=1)
    results = []
    for gene, transcripts in gene2transcripts.items():
        if len(transcripts) == 1:
            continue
        results.append([gene, pool.apply_async(pair2trans, (transcripts, transcript2exons))])
    pool.close()
    pool.join()
    #if os.path.exists("./gtf_index.txt"):
    #    return results
    return results



if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.tpm, "r") as f:
        f.readline()
        transcripts = []
        for line in f:
            #print(line)
            line = line.rstrip().split("\t")
            #transcripts.append(line[0])
            transcripts.append(line[0].split('.')[0])
        transcript2idx = dict(zip(transcripts, range(len(transcripts))))

    gtf = read_gtf(args.gtf)

    pair_results = pair_t(gtf.gene2transcripts, gtf.transcript2exons)
    #den = "exon"
    den = "site"
    with open(args.index, "w") as f:
        f.write("chr:exon_up+exon_down\tgene_id\tgene_name\texon_pari_idx\tdenominator_idx\n")
        for result in pair_results:
            gene_id = result[0]
            quant_idx2gtf_idx = dict()
            for exon_pair, isoforms in result[1].get().items():
                gene_name = gtf.geneid2name[gene_id]
                chr = gtf.gene2chr[gene_id]
                pair_idx = []
                for isoform in isoforms:
                    if isoform in transcript2idx:
                        pair_idx.append(str(transcript2idx[isoform]))
                if den == "exon":
                    den_idx = []
                    exon_0, exon_1 = exon_pair.split("+")[0].split("-"), exon_pair.split("+")[1].split("-")
                    exon_0 = (gene_id, exon_0[0], exon_0[1])
                    exon_1 = (gene_id, exon_1[0], exon_1[1])
                    transcripts_den = set(gtf.exon2transcripts[exon_0] + gtf.exon2transcripts[exon_1])
                    for transcript in transcripts_den:
                        if transcript in transcript2idx:
                            den_idx.append(str(transcript2idx[transcript]))
                if den == "site":
                    den_idx = []
                    exon_0, exon_1 = exon_pair.split("+")[0].split("-"), exon_pair.split("+")[1].split("-")
                    site_0 = (gene_id, "l", exon_0[0])
                    site_1 = (gene_id, "r", exon_0[1])
                    site_2 = (gene_id, "l", exon_1[0])
                    site_3 = (gene_id, "r", exon_1[1])
                    transcripts_den = set(gtf.site2transcripts[site_0] + gtf.site2transcripts[site_1] +
                                          gtf.site2transcripts[site_2] + gtf.site2transcripts[site_3])
                    for transcript in transcripts_den:
                        if transcript in transcript2idx:
                            den_idx.append(str(transcript2idx[transcript]))
                pair_idx.sort()
                den_idx.sort()
                if len(pair_idx) == 0 or tuple(pair_idx) == tuple(den_idx):
                    continue
                quant_idx = (tuple(pair_idx), tuple(den_idx))
                if quant_idx not in quant_idx2gtf_idx:
                    quant_idx2gtf_idx[quant_idx] = [[chr + ":" + exon_pair], gene_id, gene_name, ",".join(pair_idx), ",".join(den_idx)]
                else:
                    quant_idx2gtf_idx[quant_idx][0].append(chr+":" + exon_pair)
            for key, value in quant_idx2gtf_idx.items():
                value[0] = ";".join(value[0])
                f.write("\t".join(value) + "\n")
        f.close()
