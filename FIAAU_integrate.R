suppressMessages(library(stringr))
suppressMessages(library(getopt))
suppressMessages(library(this.path))

options<-matrix(c(
  "help", "h", "0", "logical","help",
  "prima_dir", "f", "1", "character","FIAAU dir, the output dir for FIAAU_process.py",
  "control","c","1","character","control name",
  "treatment","t","1","character","treatment name",
  "over_num","n","2","integer","method overlap number, default: 2",
  "qapa_cut","q","2","double","QAPA diff cutoff, default: 20",
  "p_met","pm","2","double","p value cutoff of each method, default: 0.05",
  "p_int","pi","2","double","p_integrate cutoff, default: 0.05"
), ncol=5, byrow=T)

args=getopt(options)

if (is.null(args$prima_dir) | is.null(args$control) | is.null(args$treatment)){
  cat(getopt(options, usage=T), "\n")
  q() 
}

if (is.null(args$over_num)) {
 args$over_num = 2
}
if (is.null(args$p_int)) {
 args$p_int = 0.05
}
if (is.null(args$p_met)) {
 args$p_met = 0.05
}
if (is.null(args$qapa_cut)) {
 args$qapa_cut = 20
}

fiaau_dir = args$prima_dir
control = args$control
treatment = args$treatment
over_num = args$over_num
p_int = args$p_int
p_met = args$p_met
qapa_cut = args$qapa_cut

check_path = function(path){
  if (grepl("^/", path)) {
    path = path  
  } else {
    path = paste0(getwd(), "/", path)	  
  }	
  return(path)
}  

print_log = function(log){
  current_time <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(paste0(current_time, " ", log, "\n"))
}

create_utr_merge_list = function(merge_utr_path){
  list1 = list()
  utr_merge = read.table(merge_utr_path,sep="\t")
  utr_merge$gene = unlist(lapply(strsplit(utr_merge$V10,"gene_name="), function(x) x[2]))
  utr_merge$gene = unlist(lapply(strsplit(utr_merge$gene,";"), function(x) x[1]))
  for(i in unique(utr_merge$gene)){
    list1[[i]] = utr_merge[utr_merge$gene == i, ]
  }
  return(list1)
}
res_type = function(res,test_method,test_thre, diff, diff_thre){
  res = res[!is.na(res[,test_method]),]
  res = res[!is.na(res[,diff]),]
  type_list = c()
  for (i in 1:nrow(res)){
    tmp = res[i,]
    if (tmp[,test_method]>=test_thre | is.na(tmp[,test_method]) | tmp[,diff] == diff_thre){
      type_list = c(type_list,"none")
    }else if(tmp[,test_method]<test_thre & tmp[,diff]>diff_thre){
      type_list = c(type_list,"long")
    }else if(tmp[,test_method]<test_thre & tmp[,diff]<(-diff_thre)){
      type_list = c(type_list,"short")
    }else{
      print(tmp)
    }
  }
  res$type = type_list
  return(res)
}
 
dapars_process = function(path, compare, p_met){
  print_log("Procssing DaPars...")
  res = read.table(paste0(path,"/DaPars/out/",compare,"_All_Prediction_Results.txt"),header = T)
  if(nrow(res) == 0){
	cat("  DaPars has no results, skip this method!\n")
	return(0)
  }
  res = res_type(res, "adjusted.P_val",p_met,"PDUI_Group_diff",0)
  res$gene = unlist(lapply(strsplit(res$Gene,"[|]"),function(x) x[2]))
  two_gene = names(table(res$gene))[table(res$gene)>1]
  flag = 1
  for(i in two_gene){
    tmp = res[res$gene==i,]
    if(nrow(tmp[tmp$adjusted.P_val<p_met,])==1){
      tmp2 = tmp[tmp$adjusted.P_val<p_met,]
    }else{
      p = min(tmp$adjusted.P_val) * nrow(tmp)
      diff = sum((1-tmp$adjusted.P_val) * tmp$PDUI_Group_diff)
      if(p < p_met){
        if(diff > 0){
          new_pred = "long"
        }else{
          new_pred = "short"
        }
      }else{
        new_pred = "none"
      }
      tmp_flt = tmp[tmp$adjusted.P_val<p_met & tmp$type == new_pred,]
      if(nrow(tmp_flt)!=0){
        strand = unlist(lapply(strsplit(tmp_flt$Gene,"[|]"), function(x) x[4]))[1]
        apa_site = paste(tmp_flt$Predicted_Proximal_APA, collapse = ",")
        tmp2 = tmp_flt[1,]
        tmp2$Predicted_Proximal_APA = apa_site
      }else{
        tmp2 = tmp[1,]
      }
      tmp2$adjusted.P_val = p
      tmp2$PDUI_Group_diff = diff
      tmp2$type = new_pred
    }
    if(flag == 1){
      df1 = tmp2
      flag = flag + 1
    }else{
      df1 = rbind(df1,tmp2)
    }
  }
  res_uniq = res[!res$gene %in% two_gene,]
  res = rbind(res_uniq,df1)
  gene_type = res[,c("gene","type","adjusted.P_val","Predicted_Proximal_APA","PDUI_Group_diff")]
  colnames(gene_type) = c("gene","pred","p","apa_site","diff")
  return(gene_type)
}
diffutr_process = function(path, compare, p_met){
  print_log("Procssing diffUTR...")
  res = read.table(paste0(path,"/diffUTR/out/",compare,".txt"),header = T)
  if(nrow(res) == 0){
	cat("  diffUTR has no results, skip this method!\n")
	return(0)
  }
  res = res_type(res, "q.value",p_met,"w.coef",0)  
  gene_type = res[,c("name","type", "q.value","w.coef")]
  colnames(gene_type) = c("gene","pred","p","diff")
  return(gene_type)
}

labrat_process = function(path, compare, ensg, p_met){
  print_log("procssing LABRAT...")
  rownames(ensg) = ensg$V2
  res = read.table(paste0(path,"/LABRAT/out/",compare,".LABRAT.psis.pval"),header = T)
  if(nrow(res) == 0){
	cat("  LABRAT has no results, skip this method!\n")
	return(0)
  }
  #res = res[res$genetype == "TUTR",]
  res = res_type(res, "FDR",p_met,"deltapsi",0)
  gene_type = res[,c("Gene","type", "FDR","deltapsi")]
  gene_list = unlist(lapply(strsplit(gene_type$Gene,"[.]"), function(x) x[1]))
  gene_list = ensg[gene_list,1]
  gene_type$Gene = gene_list
  colnames(gene_type) = c("gene","pred","p","diff")
  return(gene_type)
}

apatrap_process = function(path, compare, trans, p_met){
  print_log("procssing APAtrap...")
  trans$V1 = unlist(lapply(strsplit(trans$V1,"[.]"), function(x) x[1]))
  rownames(trans) = trans$V1
  res = read.table(paste0(path,"/APAtrap/out/",compare,"_deAPA.txt"),header = T)
  if(nrow(res) == 0){
	cat("  APAtrap has no results, skip this method!\n")
	return(0)
  }
  res = res_type(res, "p.adjust",p_met,"r",0)  
  gene_type = res[,c("Gene","type", "p.adjust","Predicted_APA","r")]
  colnames(gene_type) = c("trans","pred","p","APA","r")
  gene_type$strand = unlist(lapply(strsplit(gene_type$trans,"[|]"), function(x) x[4]))
  gene_type$trans = unlist(lapply(strsplit(gene_type$trans,"[.]"), function(x) x[1]))
  gene_type$gene = trans[gene_type$trans,2]
  flag = 1
  for(i in unique(gene_type$gene)){
    tmp_res = gene_type[gene_type$gene==i,]
    if(nrow(tmp_res)==1){
      tmp_r_df = data.frame(gene=i,r=tmp_res$r)
      tmp_p_df = data.frame(gene=i,p=tmp_res$p)
    }else{
      if(nrow(tmp_res[tmp_res$p<p_met,]) == 1){
        diff = tmp_res[tmp_res$p<p_met,"r"]
        p = tmp_res[tmp_res$p<p_met,"p"]
        tmp_r_df = data.frame(gene = i, r=diff)
        tmp_p_df = data.frame(gene = i, p=p)
      }else{
        p <- min(tmp_res$p) 
        diff = sum((1-tmp_res$p) * tmp_res$r)
        tmp_r_df = data.frame(gene = i, r=diff)
        tmp_p_df = data.frame(gene = i, p=p)
      }
    }
    if(flag == 1){
      gene_type_r = tmp_r_df
      gene_type_p = tmp_p_df
      flag = 2
    }else{
      gene_type_r = rbind(gene_type_r, tmp_r_df)
      gene_type_p = rbind(gene_type_p, tmp_p_df)
    }
  }
  rownames(gene_type_r) = gene_type_r$gene
  rownames(gene_type_p) = gene_type_p$gene
  gene_type_p = gene_type_p[gene_type_r$gene,]
  gene_type_r$p = gene_type_p$p
  res2 = res_type(gene_type_r, "p",p_met,"r",0)
  colnames(res2)[colnames(res2)=="type"] = "pred"
  colnames(res2) = c("gene","diff","p","pred")
  return(res2)
}

qapa_process = function(path, compare, merge_utr_list, qapa_cut){
  print_log("procssing QAPA...")
  res = read.table(paste0(path,"/QAPA/out/",compare,"_pau_result.txt"),header = T,sep = "\t")
  if(nrow(res) == 0){
	cat("  QAPA has no results, skip this method!\n")
	return(0)
  }
  res = res[!is.na(res$DPAU),]
  res = res[!is.na(res$Gene),]
  col_num = ncol(res)
  type_list = c()
  p_d_diff_list = c()
  gene_list = c()
  for (i in unique(res$Gene_Name)){
    tmp_utr = merge_utr_list[[i]]
    tmp = res[res$Gene_Name==i,]
    gene_list = c(gene_list, i)
    apa = unlist(lapply(strsplit(tmp[,1],"_"),function(x) x[3]))
    diff = 0
    if (length(apa)==1){
      if (apa=="S"){
        if(tmp[,"DPAU"]>qapa_cut){
          type_list = c(type_list,"long")
        }else if(tmp[,"DPAU"]< (-qapa_cut)){
          type_list = c(type_list,"short")
        }else{
          type_list = c(type_list,"none")
        }
        p_d_diff_list = c(p_d_diff_list, tmp[1,"DPAU"])
      }
    }else{
      tmp = tmp[tmp$DPAU!=0,]
      if(nrow(tmp)==0){
        type_list = c(type_list,"none")
        p_d_diff_list = c(p_d_diff_list, 0)
      }else if(nrow(tmp)==1){
        print(i)
        single_diff = tmp[,"DPAU"]
        if(apa == "D"){
          if(single_diff <(-qapa_cut)){
            type = "short"
          }else if(single_diff > qapa_cut){
            type = "long"
          }else{
            type = "none"
          }
        }else if(apa == "P"){
          if(single_diff<(-qapa_cut)){
            type = "long"
          }else if(single_diff > qapa_cut){
            type = "short"
          }else{
            type = "none"
          }
        }
        p_d_diff_list = c(p_d_diff_list, single_diff)
        type_list = c(type_list,type)
      }else{
        strand = unique(tmp$Strand)
        if(strand == "+"){
          sort_idx = order(tmp[,10])
        }else{
          sort_idx = order(tmp[,9])
        }
        tmp = tmp[sort_idx,]
        comb_mat = t(combn(nrow(tmp), 2))
        comb_num = 0
        for(j in 1:nrow(comb_mat)){
          tmp_comb = comb_mat[j,]
          if (strand == "+"){
            tmp_prox = tmp[tmp_comb[1],"DPAU"]
            tmp_dist = tmp[tmp_comb[2],"DPAU"]
            prox_site = tmp[tmp_comb[1],10]
            dist_site = tmp[tmp_comb[2],10]
          }else if(strand == "-"){
            tmp_prox = tmp[tmp_comb[2],"DPAU"]
            tmp_dist = tmp[tmp_comb[1],"DPAU"]
            prox_site = tmp[tmp_comb[2],9]
            dist_site = tmp[tmp_comb[1],9]
          }
          for(k in 1:nrow(tmp_utr)){
            tmp_utr2 = tmp_utr[k,]
            utr_s = tmp_utr2[,2]
            utr_e = tmp_utr2[,3]
            flag = 0
            if(prox_site >= utr_s & prox_site <= utr_e & dist_site >= utr_s & dist_site <= utr_e){
              comb_num = comb_num + 1
              tmp_diff = tmp_dist - tmp_prox
              diff = diff + tmp_diff
            }
          }
        }
        if (comb_num != 0){
          p_d_diff_list = c(p_d_diff_list, diff/(2*comb_num))
          if (diff/comb_num > qapa_cut){
            type_list = c(type_list,"long")
          }else if(diff/comb_num < (-qapa_cut)){
            type_list = c(type_list,"short")
          }else{
            type_list = c(type_list,"none")
          }
        }else{
          p_d_diff_list = c(p_d_diff_list, 0)
          type_list = c(type_list,"none")
        }
      }
    }
  }
  gene_type = data.frame(gene=gene_list,pred=type_list,diff=p_d_diff_list,p=p_d_diff_list)
  return(gene_type)
}

csi_process = function(path,compare, merge_utr_list, csi_anno, p_met){
  print_log("Procssing CSI_UTR...")
  csi_compare = gsub("vs","VS",compare)
  res = read.table(paste0(path,"/CSI_UTR/out/DifferentialExpression/WITHIN_UTR/",csi_compare,"_CSI.WITHINUTR.diff.txt"),header = T,sep = "\t")
  if(nrow(res) == 0){
	cat("  CSI_UTR has no results, skip this method!\n")
	return(0)
  }
  csi_anno[,"gene"] = csi_anno$V8
  gene_df = csi_anno[,c("gene","V6")]
  gene_df = unique(gene_df)
  rownames(gene_df) = gene_df$gene
  type_list = c()
  p_list = c()
  diff_list = c()
  gene_list = c()
  for (i in unique(res$GENE_SYM)){
    tmp_utr = merge_utr_list[[i]]
    gene_list = c(gene_list, i)
    tmp = res[res$GENE_SYM==i,]
    strand = gene_df[i,2]
    sort_idx = order(unlist(lapply(strsplit(tmp[,1],"[:-_]"), function(x) x[4])))
    tmp = tmp[sort_idx,]
    list1 = list()
    diff = 0
    if(nrow(tmp) == 1){
      if(tmp$FDR>=p_met){
        type_list = c(type_list,"none")
      }else{
        if(tmp[,ncol(tmp)-2] < 0){
          type_list = c(type_list,"short")
        }else if(tmp[,ncol(tmp)-2] > 0){
          type_list = c(type_list,"long")
        }else{
          type_list = c(type_list,"none")
        }
      }
      p_list = c(p_list, tmp$FDR)
      diff_list = c(diff_list, tmp[,ncol(tmp)-2])
    }else{
      p_comb <- min(tmp$FDR)
      comb_mat = t(combn(nrow(tmp), 2))
      diff = 0
      list1 = c()
      comb_num = 0
      for(j in 1:nrow(comb_mat)){
        tmp_comb = comb_mat[j,]
        if (strand == "+"){
          tmp_prox = tmp[tmp_comb[1],ncol(tmp)-2]
          tmp_dist = tmp[tmp_comb[2],ncol(tmp)-2]
          tmp_prox_p = tmp[tmp_comb[1],ncol(tmp)]
          tmp_dist_p = tmp[tmp_comb[2],ncol(tmp)]
          tmp_prox_site = strsplit(tmp[tmp_comb[1],1], "-")[[1]][2]
          tmp_dist_site = strsplit(tmp[tmp_comb[2],1], "-")[[1]][2]
        }else if(strand == "-"){
          tmp_prox = tmp[tmp_comb[2],ncol(tmp)-2]
          tmp_dist = tmp[tmp_comb[1],ncol(tmp)-2]
          tmp_prox_p = tmp[tmp_comb[2],ncol(tmp)]
          tmp_dist_p = tmp[tmp_comb[1],ncol(tmp)]
          tmp_prox_site = strsplit(tmp[tmp_comb[2],1], "-")[[1]][2]
          tmp_dist_site = strsplit(tmp[tmp_comb[1],1], "-")[[1]][2]
        }
        for(k in 1:nrow(tmp_utr)){
          tmp_utr2 = tmp_utr[k,]
          utr_s = tmp_utr2[,2]
          utr_e = tmp_utr2[,3]
          flag = 0
          if(tmp_prox_site >= utr_s & tmp_prox_site <= utr_e & tmp_dist_site >= utr_s & tmp_dist_site <= utr_e){
            comb_num = comb_num + 1
            tmp_diff = tmp_dist*(1-tmp_dist_p) - tmp_prox * (1-tmp_prox_p)
            diff = diff + tmp_diff
          }
        }
      }
      p_list = c(p_list, p_comb)
      if(comb_num != 0){
        diff_list = c(diff_list, diff/(2*comb_num))
        if(p_comb<p_met){
          if (diff > 0){
            type_list = c(type_list,"long")
          }else if(diff < 0){
            type_list = c(type_list,"short")
          }else{
            type_list = c(type_list,"none")
          }
        }else{
          type_list = c(type_list,"none")
        }
      }else{
        diff_list = c(diff_list, 0)
        type_list = c(type_list,"none")
      }
    }
  }
  gene_type = data.frame(gene=gene_list,pred=type_list,p=p_list,diff=diff_list)
  return(gene_type)
}
combp = function(p_list){
  p_num = length(p_list)
  x = -2*sum(log(p_list))
  p_comb = 1-pchisq(x,2*p_num)
  method_w = 0
  for(i in 1:p_num){
    method_w = method_w + 1/i
  }
  p_w = p_comb
  return(p_w)
}
combc = function(p_list, change_list){
  change_point = weighted.mean(change_list, w = 1 - p_list, na.rm = TRUE)
  if(change_point==0){
    change = "none"
  }else if(change_point>0){
    change = "long"
  }else if(change_point<0){
    change = "short"
  }else{
    change = "none"
  }
  res_list = list(cp=change_point, change = change)
  return(res_list)
}
combres = function(res_list, comb_method_list, p_int, over_num){
  method_list = names(res_list)
  comb_method_list = comb_method_list[comb_method_list %in% method_list]
  gene_union = c()
  for(i in comb_method_list){
    gene_union = c( gene_union, res_list[[i]][,"gene"])
    if(i == "dapars"){
      dapars_res = res_list[[i]]
    }else if(i == "apatrap"){
      apatrap_res = res_list[[i]]
    }
  }
  gene_union = unique(gene_union)
  gene_num = length(gene_union)
  flag = 1
  for(i in gene_union){
    p_list = c()
    change_list = c()
    method_list = c()
    method_type_list = c()
    method_change_list = c()
    w_list = c()
    for(j in comb_method_list){
      tmp_res = res_list[[j]]
      if(i %in% tmp_res$gene){
        gene_tmp_res = tmp_res[tmp_res$gene==i,]
        if(j != "qapa"){
          p_list = c(p_list, gene_tmp_res[,"p"])
          w_list = c(w_list, gene_tmp_res[,"p"])
        }else{
          w_list = c(w_list, abs(gene_tmp_res[,"p"])/100)
        }
        change = gene_tmp_res[,"pred"]
        if(change == "none"){
          change_list = c(change_list,0)
        }else if(change == "short"){
          change_list = c(change_list,-1)
        }else if(change == "long"){
          change_list = c(change_list,1)
        }
        method_list = c(method_list, j)
        method_type_list = c(method_type_list, paste0(j,"=",change))
        method_change_list = c(method_change_list,change)
      }
    }
    if(length(p_list)!=0){
      p_comb = combp(p_list)
      change_res = combc(w_list, change_list)
      pred = change_res$change
      if(p_comb>=p_int){
        pred = "none"
      }
    }else{
      change_res = combc(w_list, change_list)
      pred = change_res$change
      if(pred=="none"){
        p_comb = 2
      }else{
        p_comb = -1
      }
    }
    pred_over = sum(method_change_list == pred)
    if(pred_over < over_num){
      pred = "none"
    }
    tmp_df=data.frame(gene=i,pred=pred,p_int=p_comb,score=change_res$cp,method=length(method_list),
                      method_over=pred_over,method=paste(method_list,collapse = ","),
                      method_type=paste(method_type_list,collapse = ","))
    if(flag == 1){
      res_df = tmp_df
      flag = flag + 1
    }else{
      res_df = rbind(res_df, tmp_df)
      flag = flag + 1
    }
  }
  return(res_df)
}
create_num_df = function(res_list, method_input, pred_key, gene_anno){
  method_list = c()
  type_list = c()
  num_list = c()
  for (i in method_input){
    cat(paste0(i,"\n"))
    tmp = res_list[[i]]
    all_num = nrow(tmp)
    tb = table(tmp[,pred_key])
    for (j in names(tb)){
      method_list = c(method_list,i)
      type_list = c(type_list,j)
      num_list = c(num_list,tb[j])
    }
  }
  apa_num_df = data.frame(method=method_list, type=type_list, num=num_list)
  apa_num_df$method = factor(apa_num_df$method, levels = unique(apa_num_df$method))
  return(apa_num_df)
}
				   
fiaau_dir = check_path(fiaau_dir)
out_dir = paste0(fiaau_dir,"/FIAAU_integrate/")
if(!dir.exists(out_dir)){
  print_log(paste0("Create FIAAU integrate file: ", fiaau_dir, "/FIAAU_integrate"))
  dir.create(out_dir)
}

print_log("Load annotation data")
current_directory <- this.dir()
merge_utr_list = readRDS(paste0(current_directory, "/data/utr_merge.rds"))
ensg = readRDS(paste0(current_directory, "/data/ensg_name.rds"))
trans = readRDS(paste0(current_directory, "/data/trans_name.rds"))
csi_anno = read.table(paste0(current_directory, "/apps/CSI-UTR/CSI-UTR_v1.1.0/data/annotations/Mm10.CSIs.annot.bed"), sep = "\t")

comb_method_list = c("diffutr","dapars","labrat","apatrap","qapa","csiutr")
compare = paste0(treatment, "_vs_", control)
tmp_list = list(
	diffutr = diffutr_process(fiaau_dir, compare, p_met),
	dapars = dapars_process(fiaau_dir, compare, p_met),
	labrat = labrat_process(fiaau_dir, compare, ensg, p_met),
	apatrap = apatrap_process(fiaau_dir, compare, trans, p_met),
	qapa = qapa_process(fiaau_dir, compare, merge_utr_list, qapa_cut),
	csiutr=csi_process(fiaau_dir,compare,merge_utr_list, csi_anno, p_met)
)
res_list = list()
for(i in comb_method_list){
	if(!is.null(nrow(tmp_list[[i]]))){
		res_list[[i]] = tmp_list[[i]]
	}else{
		comb_method_list = comb_method_list[comb_method_list != i]
	}
}

print_log(paste0("Integrate ", length(comb_method_list), " (", paste(sort(comb_method_list), collapse=","), ") results..."))
res_list$comb = combres(res_list, comb_method_list, p_int, over_num)
write.table(res_list$comb,paste0(out_dir,"/",compare,"_integrate_result.txt"),sep="\t", col.names = T, row.names = F,quote = F)
print_log("Finished")
