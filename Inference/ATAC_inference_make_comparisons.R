coldata <- readRDS("coldata.rds")

coldata <- coldata[!(is.na(coldata$purity) | coldata$purity == 0 | coldata$subtissue == "adenoma"), ]
coldata$group <- NA

#525: all region C vs rest
coldata[coldata$patient=="C525", "group"] <- "C525.ABD"
coldata[coldata$patient=="C525" & coldata$region =="C", "group"] <- "C525.C"
#538: all region D vs rest
coldata[coldata$patient=="C538" & coldata$region !="A" , "group"] <- "C538.BC"
coldata[coldata$patient=="C538" & coldata$region =="D", "group"] <- "C538.D"
#539: 3 groups: allA+B1_G3 vs C+D vs rest
coldata[coldata$patient=="C539", "group"] <- "C539.B"
coldata[coldata$patient=="C539" & coldata$region =="A" | row.names(coldata)=="C539_B1_G3", "group"] <- "C539.AB1_G3"
coldata[coldata$patient=="C539" & coldata$region =="C" | coldata$patient=="C539" & coldata$region =="D" , "group"] <- "C539.CD"
#542: 17/05 ABD vs. C
coldata[coldata$patient=="C542" & coldata$region !="F", "group"] <- "C542.ABD"
coldata[coldata$patient=="C542" & coldata$region=="C", "group"] <- "C542.C"
#549: region A vs rest
coldata[coldata$patient=="C549" & coldata$region =="A", "group"] <- "C549.A"
coldata[coldata$patient=="C549" & coldata$region !="A", "group"] <- "C549.BCD"
#516 B vs A
coldata[coldata$patient=="C516" & coldata$region =="B", "group"] <- "C516.B"
coldata[coldata$patient=="C516" & coldata$region =="A", "group"] <- "C516.A"
#518 D vs A+B vs C
coldata[coldata$patient=="C518" & coldata$region =="D", "group"] <- "C518.D"
coldata[coldata$patient=="C518" & coldata$region =="C", "group"] <- "C518.C"
coldata[coldata$patient=="C518" & coldata$region =="A", "group"] <- "C518.AB"
coldata[coldata$patient=="C518" & coldata$region =="B", "group"] <- "C518.AB"
#524 B vs rest
coldata[coldata$patient=="C524", "group"] <- "C524.ACD"
coldata[coldata$patient=="C524" & coldata$region =="B", "group"] <- "C524.B"
#531 B vs rest
coldata[coldata$patient=="C531", "group"] <- "C531.ACD"
coldata[coldata$patient=="C531" & coldata$region =="B", "group"] <- "C531.B"
#551 all vs C1_G4+B1_G3+A1_G6+B1_G7+B1_G2+A1_G9
coldata[coldata$patient=="C551", "group"] <- "C551.rest"
coldata[row.names(coldata) %in% paste0("C551_",c("C1_G4","B1_G3","A1_G6","B1_G7","B1_G2","A1_G9")), "group"] <- "C551.select"
#559 B vs rest
coldata[coldata$patient=="C559", "group"] <- "C559.CD"
coldata[coldata$patient=="C559" & coldata$region =="B", "group"] <- "C559.B"
#562 all vs A1_G7
coldata[coldata$patient=="C562", "group"] <- "C562.ABCD"
coldata[row.names(coldata)=="C562_A1_G7", "group"] <- "C562.A1_G7"

saveRDS(coldata,"coldata_inf.rds")
