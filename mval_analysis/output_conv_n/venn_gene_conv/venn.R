library("VennDiagram")

install.packages("VennDiagram")

draw.pairwise.venn(area1 = 73,
                   area2 = 16,
                   cross.area = 15,
                   category = c("ctrl_worse" , "em_worse"))


draw.quintuple.venn(
  area1 = 60,
  area2 = 118,
  area3 = 111,
  area4 = 124,
  area5 = 90,
  n12 = 1, n13 = 3, n14 = 1, n15 = 3,
  n23 = 8, n24 = 7, n25 = 10,
  n34 = 9, n35 = 3, n45 = 4,
  n123 = 0, n124 = 0, n125 = 0,
  n134 = 0, n135 = 0, n145 = 0,
  n234 = 0, n235 = 0, n245 = 2,
  n345 = 1,
  n1234 = 0, n1235 = 0, n1245 = 0, n1345 = 0, n2345 = 0,
  n12345 = 0,
  category = c("exp.G1" , "exp.G2", "exp.G3", "exp.G4", "exp.G5"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"))
