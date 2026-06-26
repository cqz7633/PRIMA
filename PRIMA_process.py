#coding=utf8
import os
import argparse
import sys
from datetime import datetime

from scripts.utils import *

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="FIAAU for 6 methods processed")
parser.add_argument("-f", type=str, default="", help="absolute paths of fastq information file")
parser.add_argument("-b", type=str, default="", help="absolute paths of bam information file")
parser.add_argument("-c", type=str, default="", help="control")
parser.add_argument("-t", type=str, default="", help="treatment")
parser.add_argument("-p", type=str, default="y", help="paired or single ('y' or 'n') default: y")
parser.add_argument("-m", type=int, default="4", help="core numbers default: 4")
parser.add_argument("-r", type=int, default="150", help="reads length  default: 150")
parser.add_argument("-o", type=str, default="/PRIMA", help="out put dir  default: ./PRIMA_Y-m-d_H-M-S")
parser.add_argument("-l", type=str, default="n", help="run the step 'makeTFfasta' of LABRAT or not ('y' or 'n') default: n")
parser.add_argument("-bs", type=int, default="10", help="bin size for calculating big wig file default: 10")
parser.add_argument("-ct", type=float, default="0.5", help="coverage cut off default: 0.5")

args = parser.parse_args()

# get current time
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d_%H-%M-%S")
current_dir = os.getcwd()

# args var
fq_table = args.f
bam_table = args.b
paired = args.p
core = args.m
reads_len = args.r
control = args.c
treatment = args.t
maketf = args.l
bin_size = args.bs
cov_cutoff = args.ct
if args.o == "/FIAAU":
    out_path = current_dir + args.o + "_" + formatted_time + "/"
else:
    out_path = check_path(args.o)

# args judge
if fq_table == "":
    fq_war1 = "fastq table is not input!"
    sys.exit(fq_war1)
if bam_table == "":
    bam_war1 = "bam table is not input!"
    sys.exit(bam_war1)
if reads_len == "":
    rlen_war1 = "reads length is not input!"
    sys.exit(rlen_war1)
if not isinstance(reads_len, int):
    rlen_war2 = "reads length (-r) is not correct!\nreads length is must be int number!"
    sys.exit(rlen_war2)
if control == "" or treatment == "":
    ct_war = "control/treatment is not input!"
    sys.exit(ct_war)
if not isinstance(core, int):
    core_war = "core (-c) is not correct!\ncore is must be int number!"
    sys.exit(core_war)
if not os.path.exists(bam_table):
    btb_war = "bam_table file is not exists!"
    sys.exit(btb_war)
if not os.path.exists(fq_table):
    ftb_war = "fq_table file is not exists!"
    sys.exit(ftb_war)
if not os.path.exists(out_path):
    print(out_path + " is not exists, creating ...")
    os.makedirs(out_path)
if not isinstance(bin_size, int):
    rlen_war2 = "reads length (-r) is not correct!\nreads length is must be int number!"
    sys.exit(rlen_war2)
if not isinstance(cov_cutoff, float):
    rlen_war2 = "coverage cut off (-ct) is not correct!\ncoverage cut off is must be float type!"
    sys.exit(rlen_war2)
    
script_dir = os.path.dirname(os.path.abspath(__file__))
#====================== annotation =========================#
#dapars anno
dapars_bed = script_dir + "/anno/mm10/dapars_anno.bed"
csiutr_bed = script_dir + "/anno/mm10/csi.bed"
csiutr_anno = script_dir + "/anno/mm10/csi_anno.bed"
csiutr_anno_dir = script_dir + "/apps/CSI-UTR/CSI-UTR_v1.1.0/data"
#qapa anno
qapa_genome_fa = script_dir + "/anno/mm10/mm10.fa"
qapa_utr = script_dir + "/anno/mm10/qapa_anno.bed"
qapa_ident = script_dir + "/anno/mm10/qapa_ident.txt"
#apatrap anno
apatrap_genemodel = script_dir + "/anno/mm10/apatrap_anno.bed"
apatrap_bed = script_dir + "/anno/mm10/apatrap_3utr.bed"
#diffutr anno
diffutr_bed = script_dir + "/anno/mm10/diffutr_anno.bed"
#labrat anno
labrat_gff = script_dir + "/anno/mm10/labrat_anno.addhead.gff3"
labrat_fa = script_dir + "/anno/mm10/mm10.fa"
labrat_seq = script_dir + "/anno/mm10/TFseqs.fasta"
labrat_db = script_dir + "/anno/mm10/labrat_anno.addhead.gff3.db"
#genomeic anno
geno_size = script_dir + "/anno/mm10/mm10.chrom.sizes"

#====================== apps =========================#
dapars_main = script_dir + "/apps/DaPars/src/DaPars_main.py"
csi_utr = script_dir + "/apps/CSI-UTR/CSI-UTR_v1.1.0/bin/CSI-UTR"
apatrap_identify = script_dir + "/apps/APAtrap/identifyDistal3UTR"
apatrap_predic = script_dir + "/apps/APAtrap/predictAPA"
apatrap_deapa = script_dir + "/scripts/run_deAPA.R"
diffutr_script = script_dir + "/scripts/run_diffutr.R"
qapa_diff_script = script_dir + "/scripts/run_qapa_diff.R"

#====================== args =========================#
csiutr_arg = "-genome Mm10 -r %s -coverage_cut %d" % (reads_len, cov_cutoff)


def main():
    bam_dict = table_cond(bam_table)
    fq_dict = table_cond(fq_table)
    if (control not in bam_dict) or (treatment not in bam_dict):
        bam_ware = "control/treatment not in bam table"
        sys.exit(bam_ware)
    if (control not in fq_dict) or (treatment not in fq_dict):
        fq_ware = "control/treatment not in fq table"
        sys.exit(fq_ware)
    cond_list = [control, treatment]
    sub_dapars_sh, dapars_res = dapars_process(dapars_main, cond_list, bam_dict, out_path, core, dapars_bed, cov_cutoff, bin_size) #in groupA long:dpdui>0
    sub_csiutr_sh, csiutr_res = csiutr_process(csi_utr, cond_list, bam_dict, out_path, csiutr_arg, csiutr_bed, csiutr_anno, csiutr_anno_dir) #long:dpsi>0
    sub_diffutr_sh, diffutr_res = diffutr_process(diffutr_script, cond_list, bam_table, out_path, diffutr_bed, "diffUTR", cov_cutoff, reads_len)
    sub_qapa_sh, qapa_res = qapa_process(fq_table, cond_list, out_path, paired, qapa_genome_fa, qapa_utr, core, qapa_ident, qapa_diff_script)
    sub_apatrap_sh, apatrap_res = apatrap_process(cond_list, bam_dict, out_path, apatrap_identify, apatrap_predic, apatrap_deapa, geno_size, apatrap_genemodel, cov_cutoff)
    sub_labrat_sh, labrat_res = labrat_process(fq_table, cond_list, out_path, labrat_fa, labrat_gff, labrat_db, labrat_seq, paired, maketf, core)
    os.system(sub_dapars_sh)
    os.system(sub_csiutr_sh)
    os.system(sub_diffutr_sh)
    os.system(sub_qapa_sh)
    os.system(sub_apatrap_sh)
    os.system(sub_labrat_sh)

if __name__ == '__main__':
    main()
