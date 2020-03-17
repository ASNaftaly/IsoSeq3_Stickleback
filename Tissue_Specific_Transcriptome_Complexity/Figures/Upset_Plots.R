#making upset plots to show overlap between single tissues and joined tissues
#Used UpSetR package in R
#Author:Alice Naftaly, March 2020

library(UpSetR)

#shared genes

allcomparisons <- c(Female.Liver=1684, Male.Liver=1323, Female.Brain=1282, Male.Brain=2244, Female.Pronephros=5120, Male.Pronephros=5251, Ovary=2017,Testis=3773, All.Female.Tissues = 7299, All.Male.Tissues = 8223, All.Somatic.Female.Tissues = 6643, All.Somatic.Male.Tissues = 7177, `Female.Liver&Male.Liver`=719, `Female.Brain&Male.Brain`=771, `Female.Pronephros&Male.Pronephros`=1679, `Ovary&Testis`=750, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros&Ovary&Testis`=203, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`=258, `Female.Liver&Female.Brain&Female.Pronephros&Ovary`=319, `Male.Liver&Male.Brain&Male.Pronephros&Testis`=417, `Female.Liver&Male.Liver&Female.Brain&Male.Brain` = 273, `Female.Liver&Male.Liver&Female.Pronephros&Male.Pronephros`=523, `Female.Liver&Male.Liver&Ovary&Testis`=330, `Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`=481, `Female.Brain&Male.Brain&Ovary&Testis`=364, `Female.Pronephros&Male.Pronephros&Ovary&Testis`=1194, `All.Female.Tissues&All.Male.Tissues`= 4375, `All.Somatic.Female.Tissues&All.Somatic.Male.Tissues`= 4449)

reduced_comparisons <- c(Female.Liver=1684, Male.Liver=1323, Female.Brain=1282, Male.Brain=2244, Female.Pronephros=5120, Male.Pronephros=5251, Ovary=2017, Testis=3773, `Female.Liver&Male.Liver`=719, `Female.Brain&Male.Brain`=771, `Female.Pronephros&Male.Pronephros`=1679, `Ovary&Testis`=750, `Female.Liver&Female.Brain&Female.Pronephros&Ovary`=319, `Male.Liver&Male.Brain&Male.Pronephros&Testis`=417, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`= 4449, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros&Ovary&Testis`=4375)

#upset plot

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedGenes_UpsetPlot_SingleTissues_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"),intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain"), list("Female.Liver","Male.Liver","Female.Pronephros", "Male.Pronephros"), list("Female.Liver","Male.Liver","Ovary","Testis"), list("Female.Brain", "Male.Brain", "Female.Pronephros","Male.Pronephros"), list("Female.Brain","Male.Brain","Ovary","Testis"), list("Female.Pronephros","Male.Pronephros","Ovary","Testis"), list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"), list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")),keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Genes",sets.x.label="Number of Genes", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedGenes_UpsetPlot_SingleTissues_Reduced_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"), intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"), list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")),keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Genes",sets.x.label="Number of Genes", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()


jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedGenes_UpsetPlot_JoinedTissues_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons),sets=c("All.Female.Tissues","All.Male.Tissues","All.Somatic.Female.Tissues","All.Somatic.Male.Tissues"),intersections=list(list("All.Somatic.Female.Tissues","All.Somatic.Male.Tissues"),list("All.Female.Tissues","All.Male.Tissues")) ,keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Genes",sets.x.label="Number of Genes", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedGenes_UpsetPlot_MainPaper_Fig_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(reduced_comparisons),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"),intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"),list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"), list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")) ,keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Genes",sets.x.label="Number of Genes", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()

#isoforms


allcomparisons_isoforms <- c(Female.Liver=2139, Male.Liver=1600, Female.Brain=1355, Male.Brain=2482, Female.Pronephros=8177, Male.Pronephros=8679, Ovary=2216,Testis=5174, All.Female.Tissues = 11653, All.Male.Tissues = 14068, All.Somatic.Female.Tissues = 10571, All.Somatic.Male.Tissues = 11432, `Female.Liver&Male.Liver`=690, `Female.Brain&Male.Brain`=719, `Female.Pronephros&Male.Pronephros`=4758, `Ovary&Testis`=1314, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros&Ovary&Testis`=158, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`=206, `Female.Liver&Female.Brain&Female.Pronephros&Ovary`=249, `Male.Liver&Male.Brain&Male.Pronephros&Testis`=320, `Female.Liver&Male.Liver&Female.Brain&Male.Brain` = 224, `Female.Liver&Male.Liver&Female.Pronephros&Male.Pronephros`=417, `Female.Liver&Male.Liver&Ovary&Testis`=253, `Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`=388, `Female.Brain&Male.Brain&Ovary&Testis`=290, `Female.Pronephros&Male.Pronephros&Ovary&Testis`=965, `All.Female.Tissues&All.Male.Tissues`= 5631, `All.Somatic.Female.Tissues&All.Somatic.Male.Tissues`= 6494)

reduced_comparisons_isoforms <- c(Female.Liver=2139, Male.Liver=1600, Female.Brain=1355, Male.Brain=2482, Female.Pronephros=8177, Male.Pronephros=8679, Ovary=2216,Testis=5174, `Female.Liver&Male.Liver`=690, `Female.Brain&Male.Brain`=719, `Female.Pronephros&Male.Pronephros`=4758, `Ovary&Testis`=1314, `Female.Liver&Female.Brain&Female.Pronephros&Ovary`=249, `Male.Liver&Male.Brain&Male.Pronephros&Testis`=320, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros`= 6494, `Female.Liver&Male.Liver&Female.Brain&Male.Brain&Female.Pronephros&Male.Pronephros&Ovary&Testis`=5631)




#upset plot

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedIsoforms_UpsetPlot_SingleTissues_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons_isoforms),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"),intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain"), list("Female.Liver","Male.Liver","Female.Pronephros", "Male.Pronephros"), list("Female.Liver","Male.Liver","Ovary","Testis"), list("Female.Brain", "Male.Brain", "Female.Pronephros","Male.Pronephros"), list("Female.Brain","Male.Brain","Ovary","Testis"), list("Female.Pronephros","Male.Pronephros","Ovary","Testis"), list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"), list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")),keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Isoforms",sets.x.label="Number of Isoforms", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedIsoforms_UpsetPlot_SingleTissues_Reduced_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons_isoforms),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"), intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"), list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")),keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Isofrms",sets.x.label="Number of Isoforms", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()


jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedIsoforms_UpsetPlot_JoinedTissues_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(allcomparisons_isoforms),sets=c("All.Female.Tissues","All.Male.Tissues","All.Somatic.Female.Tissues","All.Somatic.Male.Tissues"),intersections=list(list("All.Somatic.Female.Tissues","All.Somatic.Male.Tissues"),list("All.Female.Tissues","All.Male.Tissues")) ,keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Isoforms",sets.x.label="Number of Isoforms", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()

jpeg(file="/Users/Alice/Documents/UGA/White_Lab/Project/Iso_seq/Figures/Version1/SharedIsoforms_UpsetPlot_MainPaper_Fig_V1.jpeg",width=1600,height=1000,quality=100)

upset(fromExpression(reduced_comparisons_isoforms),sets=c("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros","Ovary", "Testis"),intersections=list(list("Female.Liver","Male.Liver"),list("Female.Brain","Male.Brain"),list("Female.Pronephros","Male.Pronephros"),list("Ovary","Testis"),list("Female.Liver","Female.Brain","Female.Pronephros","Ovary"),list("Male.Liver","Male.Brain","Male.Pronephros","Testis"),list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros"), list("Female.Liver","Male.Liver","Female.Brain","Male.Brain","Female.Pronephros","Male.Pronephros", "Ovary","Testis")) ,keep.order=TRUE,point.size=7,mainbar.y.label="Number of Shared Isoforms",sets.x.label="Number of Isoforms", line.size=1.5, text.scale=c(3.5,2.75,2.5,2,2.5,3),order.by="degree", decreasing=FALSE)

dev.off()
