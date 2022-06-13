#############################################
# list of covariates to use in modelling
myvars <- list(
  AWT = c("bc04",  "bc05",  "bc06",  "bc12",  "bc15",  "slope", "topo", "tri"),
  CAN = c("alt", "asp2", "ontprec", "ontslp", "onttemp", "ontveg", "watdist"), 
  NSW = c("cti", "disturb", "mi", "rainann", "raindq", "rugged", "soildepth", "soilfert", "solrad", "tempann", "topo"), # "vegsys"), 
  NZ = c("age", "deficit", "hillshade", "mas", "mat", "r2pet", "slope", "sseas", "toxicats", "tseas", "vpd"), 
  SA = c("sabio12", "sabio15", "sabio17", "sabio18", "sabio2", "sabio4", "sabio5", "sabio6"), 
  SWI = c("bcc", "calc", "ccc", "ddeg", "nutri", "pday", "precyy", "sfroyy", "slope", "sradyy", "swb", "topo")
)

# builing gam formulas
awt_fm <- paste("occ ~", paste(paste0("s(", myvars[["AWT"]], ")"), collapse = " + "))
can_fm <- paste("occ ~", paste(paste0("s(", myvars[["CAN"]][-6], ")"), collapse = " + "))
# can_fm <- paste(can_fm, " + ", myvars[["CAN"]][6])
nsw_fm <- occ ~ s(cti) + s(mi) + s(rainann) + s(raindq) + s(rugged) + 
  s(soildepth) + s(solrad) + s(tempann) + s(topo) + soilfert + disturb #+ vegsys
nz_fm <- occ ~ s(deficit) + s(hillshade) + s(mas) + s(mat) + s(r2pet) + 
  s(slope) + s(sseas) + s(tseas) + s(vpd)# + age
sa_fm <- paste("occ ~", paste(paste0("s(", myvars[["SA"]], ")"), collapse = " + "))
swi_fm <- paste("occ ~", paste(paste0("s(", myvars[["SWI"]][-2], ")"), collapse = " + "))
swi_fm <- paste(swi_fm, " + ", myvars[["SWI"]][2])
swi_fm <- "occ ~ s(bcc) + s(ccc) + s(ddeg) + s(nutri) + s(pday) + s(precyy) + s(slope) + s(sradyy) + s(swb) + s(topo) + sfroyy + calc"

myform <- list(AWT = awt_fm, 
               CAN = can_fm, 
               NSW = nsw_fm, 
               NZ = nz_fm, 
               SA = sa_fm, 
               SWI = swi_fm)


#############################################
# building model scopes for GLM mdeol selection
c("bc04",  "bc05",  "bc06",  "bc12",  "bc15",  "slope", "topo", "tri")
awt_scope <- list("bc04" = ~1 + bc04 + poly(bc04, 2, raw=TRUE),
                  "bc05" = ~1 + bc05 + poly(bc05, 2, raw=TRUE),
                  "bc06" = ~1 + bc06 + poly(bc06, 2, raw=TRUE),
                  "bc12" = ~1 + bc12 + poly(bc12, 2, raw=TRUE), 
                  "bc15" = ~1 + bc15 + poly(bc15, 2, raw=TRUE), 
                  "slope" = ~1 + slope + poly(slope, 2, raw=TRUE),
                  "topo" = ~1 + topo + poly(topo, 2, raw=TRUE),
                  "tri" = ~1 + tri + poly(tri, 2, raw=TRUE))

c("alt", "asp2", "ontprec", "ontslp", "onttemp", "ontveg", "watdist")
can_scope <- list("alt" = ~1 + alt + poly(alt, 2, raw=TRUE),
                  "asp2" = ~1 + asp2 + poly(asp2, 2, raw=TRUE),
                  "ontprec" = ~1 + ontprec + poly(ontprec, 2, raw=TRUE),
                  "ontslp" = ~1 + ontslp + poly(ontslp, 2, raw=TRUE), 
                  "onttemp" = ~1 + onttemp + poly(onttemp, 2, raw=TRUE),
                  "watdist" = ~1 + watdist + poly(watdist, 2, raw=TRUE),
                  "ontveg" = ~1 + ontveg)

c("cti", "disturb", "mi", "rainann", "raindq", "rugged", "soildepth", "soilfert", "solrad", "tempann", "topo", "vegsys")
nsw_scope <- list("cti" = ~1 + cti + poly(cti, 2, raw=TRUE),
                  "disturb" = ~1 + disturb + poly(disturb, 2, raw=TRUE),
                  "mi" = ~1 + mi + poly(mi, 2, raw=TRUE),
                  "rainann" = ~1 + rainann + poly(rainann, 2, raw=TRUE), 
                  "raindq" = ~1 + raindq + poly(raindq, 2, raw=TRUE), 
                  "rugged" = ~1 + rugged + poly(rugged, 2, raw=TRUE),
                  "soildepth" = ~1 + soildepth + poly(soildepth, 2, raw=TRUE),
                  "soilfert" = ~1 + soilfert + poly(soilfert, 2, raw=TRUE),
                  "solrad" = ~1 + solrad + poly(solrad, 2, raw=TRUE),
                  "tempann" = ~1 + tempann + poly(tempann, 2, raw=TRUE),
                  "topo" = ~1 + topo + poly(topo, 2, raw=TRUE))
                  # "vegsys" = ~1 + vegsys)

c("age", "deficit", "hillshade", "mas", "mat", "r2pet", "slope", "sseas", "toxicats", "tseas", "vpd")
nz_scope <- list("age" = ~1 + age,
                 "deficit" = ~1 + deficit + poly(deficit, 2, raw=TRUE),
                 "hillshade" = ~1 + hillshade + poly(hillshade, 2, raw=TRUE),
                 "mas" = ~1 + mas + poly(mas, 2, raw=TRUE), 
                 "mat" = ~1 + mat + poly(mat, 2, raw=TRUE), 
                 "r2pet" = ~1 + r2pet + poly(r2pet, 2, raw=TRUE),
                 "slope" = ~1 + slope + poly(slope, 2, raw=TRUE),
                 "sseas" = ~1 + sseas + poly(sseas, 2, raw=TRUE),
                 "tseas" = ~1 + tseas + poly(tseas, 2, raw=TRUE),
                 "vpd" = ~1 + vpd + poly(vpd, 2, raw=TRUE),
                 "toxicats" = ~1 + toxicats)

c("sabio12", "sabio15", "sabio17", "sabio18", "sabio2", "sabio4", "sabio5", "sabio6")
sa_scope <- list("sabio12" = ~1 + sabio12 + poly(sabio12, 2, raw=TRUE),
                 "sabio15" = ~1 + sabio15 + poly(sabio15, 2, raw=TRUE),
                 "sabio17" = ~1 + sabio17 + poly(sabio17, 2, raw=TRUE), 
                 "sabio18" = ~1 + sabio18 + poly(sabio18, 2, raw=TRUE), 
                 "sabio2" = ~1 + sabio2 + poly(sabio2, 2, raw=TRUE),
                 "sabio4" = ~1 + sabio4 + poly(sabio4, 2, raw=TRUE),
                 "sabio5" = ~1 + sabio5 + poly(sabio5, 2, raw=TRUE),
                 "sabio6" = ~1 + sabio6 + poly(sabio6, 2, raw=TRUE))

c("bcc", "calc", "ccc", "ddeg", "nutri", "pday", "precyy", "sfroyy", "slope", "sradyy", "swb", "topo")
swi_scope <- list("bcc" = ~1 + bcc + poly(bcc, 2, raw=TRUE),
                  "ccc" = ~1 + ccc + poly(ccc, 2, raw=TRUE),
                  "ddeg" = ~1 + ddeg + poly(ddeg, 2, raw=TRUE), 
                  "nutri" = ~1 + nutri + poly(nutri, 2, raw=TRUE), 
                  "pday" = ~1 + pday + poly(pday, 2, raw=TRUE),
                  "precyy" = ~1 + precyy + poly(precyy, 2, raw=TRUE),
                  "sfroyy" = ~1 + sfroyy + poly(sfroyy, 2, raw=TRUE),
                  "slope" = ~1 + slope + poly(slope, 2, raw=TRUE),
                  "sradyy" = ~1 + sradyy + poly(sradyy, 2, raw=TRUE),
                  "swb" = ~1 + swb + poly(swb, 2, raw=TRUE),
                  "topo" = ~1 + topo + poly(topo, 2, raw=TRUE),
                  "calc" = ~1 + calc)


myscope <- list(AWT = awt_scope,
                CAN = can_scope,
                NSW = nsw_scope,
                NZ = nz_scope,
                SA = sa_scope,
                SWI = swi_scope)

#############################################