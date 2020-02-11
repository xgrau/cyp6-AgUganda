library(pheatmap)
library(gmodels)


# input

tau = read.table("results_allpop/cyp6_all_cyp6.allele_fq.csv", sep="\t", header = T)


rownames(tau) = paste(tau$gene_eff, tau$PEP_eff, tau$POS)
varbool = tau$all>0.05 & tau$is_nonsyn == "True"
tam = tau[varbool,c("AOM","BFM","BFS","CMS","GAS","GWA","KES","UGS")]

col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
cod.fun = colorRampPalette(c("firebrick4","orangered", "floralwhite", "deepskyblue","dodgerblue4"))

pdf(file=paste("results_allpop/cyp6_all_cyp6.allele_fq_large.pdf",sep=""),height=30,width=12)
pheatmap(tam, color = col.fun(20), breaks = seq(0,1,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_color = "red",
         main=paste("nonsyn freqs"))
dev.off()



