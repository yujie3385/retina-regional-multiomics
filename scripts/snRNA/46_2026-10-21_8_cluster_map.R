# 2026-10-21_8 整合对象(res 0.4, 36 clusters) -> celltype 映射
# 依据: FindAllMarkers top10 + 经典marker平均表达 + 原注释构成(人工核对)
# 边界情况: cl6(VSX2 1.97/CABP5 2.90, AC标记≈0)=BC亚型(SCG2/BHLHE22程序);
#           cl29(TFAP2B/SLC6A9/CALB1)=AC; cl30(RHO 168)=Rod; cl35(ARR3/RCVRN高)=Cone
cluster_map <- c(
  "0" = "Rod",  "1" = "Rod",  "2" = "Rod",  "12" = "Rod", "17" = "Rod", "30" = "Rod",
  "10" = "Cone", "35" = "Cone",
  "4" = "BC",  "6" = "BC",  "8" = "BC",  "9" = "BC",  "15" = "BC", "16" = "BC",
  "19" = "BC", "20" = "BC", "22" = "BC", "23" = "BC", "24" = "BC", "27" = "BC", "28" = "BC",
  "18" = "HC",
  "13" = "AC", "21" = "AC", "29" = "AC", "31" = "AC", "32" = "AC", "34" = "AC",
  "25" = "RGC",
  "3" = "MG",  "5" = "MG",  "7" = "MG",  "11" = "MG", "14" = "MG", "26" = "MG", "33" = "MG"
)
