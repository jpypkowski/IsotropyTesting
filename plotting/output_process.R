# SIMULATION OUTPUT PROCESSING
# Getting p-values and rejection rates

setwd("IsotropyTesting")


pval <- function(x) (sum(as.numeric(x[1]) < x[2:1001], na.rm=TRUE)+1)/(1+sum(!is.na(x[2:1001])))
pval_SR <- function(x) (sum(as.numeric(x[1]) < x[2:100], na.rm=TRUE)+1)/(1+sum(!is.na(x[2:100])))

#########
# LGCP
LLS1 <- read.csv("LGCP\\LGCP_LGCP_1_SMALL.csv")[,-1]
LLS08 <- read.csv("LGCP\\LGCP_LGCP_08_SMALL.csv")[,-1]
LLS06 <- read.csv("LGCP\\LGCP_LGCP_06_SMALL.csv")[,-1]
LLS04 <- read.csv("LGCP\\LGCP_LGCP_04_SMALL.csv")[,-1]


LLS1_p <- rbind(apply(LLS1[2:1002,], 2, pval),
                apply(LLS1[2006:3006,], 2, pval),
                apply(LLS1[4010:5010,], 2, pval))
LLS08_p <- rbind(apply(LLS08[2:1002,], 2, pval),
                apply(LLS08[2006:3006,], 2, pval),
                apply(LLS08[4010:5010,], 2, pval))
LLS06_p <- rbind(apply(LLS06[2:1002,], 2, pval),
                apply(LLS06[2006:3006,], 2, pval),
                apply(LLS06[4010:5010,], 2, pval))
LLS04_p <- rbind(apply(LLS04[2:1002,], 2, pval),
                apply(LLS04[2006:3006,], 2, pval),
                apply(LLS04[4010:5010,], 2, pval))
POW_LLS <- cbind(apply(LLS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLS04_p, 1, function(x) {mean(x<=0.05)}))

LTS1 <- read.csv("LGCP\\LGCP_THOMAS_1_SMALL_long_short.csv")[,-1]
LTS08 <- read.csv("LGCP\\LGCP_THOMAS_08_SMALL_long_short.csv")[,-1]
LTS06 <- read.csv("LGCP\\LGCP_THOMAS_06_SMALL_long_short.csv")[,-1]
LTS04 <- read.csv("LGCP\\LGCP_THOMAS_04_SMALL_long_short.csv")[,-1]

LTS1_p <- rbind(apply(LTS1[2:1002,], 2, pval),
                apply(LTS1[2006:3006,], 2, pval),
                apply(LTS1[4010:5010,], 2, pval))
LTS08_p <- rbind(apply(LTS08[2:1002,], 2, pval),
                 apply(LTS08[2006:3006,], 2, pval),
                 apply(LTS08[4010:5010,], 2, pval))
LTS06_p <- rbind(apply(LTS06[2:1002,], 2, pval),
                 apply(LTS06[2006:3006,], 2, pval),
                 apply(LTS06[4010:5010,], 2, pval))
LTS04_p <- rbind(apply(LTS04[2:1002,], 2, pval),
                 apply(LTS04[2006:3006,], 2, pval),
                 apply(LTS04[4010:5010,], 2, pval))
POW_LTS <- cbind(apply(LTS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTS04_p, 1, function(x) {mean(x<=0.05)}))

LSS1 <- read.csv("LGCP\\LGCP_SR_1_SMALL.csv")[,-1]
LSS08 <- read.csv("LGCP\\LGCP_SR_08_SMALL.csv")[,-1]
LSS06 <- read.csv("LGCP\\LGCP_SR_06_SMALL.csv")[,-1]
LSS04 <- read.csv("LGCP\\LGCP_SR_04_SMALL.csv")[,-1]

LSS1_p <- rbind(apply(LSS1[2:101,], 2, pval_SR),
                apply(LSS1[204:305,], 2, pval_SR),
                apply(LSS1[406:505,], 2, pval_SR))
LSS08_p <- rbind(apply(LSS08[2:101,], 2, pval_SR),
                 apply(LSS08[204:305,], 2, pval_SR),
                 apply(LSS08[406:505,], 2, pval_SR))
LSS06_p <- rbind(apply(LSS06[2:101,], 2, pval_SR),
                 apply(LSS06[204:305,], 2, pval_SR),
                 apply(LSS06[406:505,], 2, pval_SR))
LSS04_p <- rbind(apply(LSS04[2:101,], 2, pval_SR),
                 apply(LSS04[204:305,], 2, pval_SR),
                 apply(LSS04[406:505,], 2, pval_SR))
POW_LSS <- cbind(apply(LSS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSS04_p, 1, function(x) {mean(x<=0.05)}))

LT5S1 <- read.csv("LGCP\\LGCP_1_TILING_5_SMALL.csv")[,-1]
LT5S08 <- read.csv("LGCP\\LGCP_08_TILING_5_SMALL.csv")[,-1]
LT5S06 <- read.csv("LGCP\\LGCP_06_TILING_5_SMALL.csv")[,-1]
LT5S04 <- read.csv("LGCP\\LGCP_04_TILING_5_SMALL.csv")[,-1]

LT5S1_p <- rbind(apply(LT5S1[2:1002,], 2, pval),
                apply(LT5S1[2006:3006,], 2, pval),
                apply(LT5S1[4010:5010,], 2, pval))
LT5S08_p <- rbind(apply(LT5S08[2:1002,], 2, pval),
                 apply(LT5S08[2006:3006,], 2, pval),
                 apply(LT5S08[4010:5010,], 2, pval))
LT5S06_p <- rbind(apply(LT5S06[2:1002,], 2, pval),
                 apply(LT5S06[2006:3006,], 2, pval),
                 apply(LT5S06[4010:5010,], 2, pval))
LT5S04_p <- rbind(apply(LT5S04[2:1002,], 2, pval),
                 apply(LT5S04[2006:3006,], 2, pval),
                 apply(LT5S04[4010:5010,], 2, pval))
POW_LT5S <- cbind(apply(LT5S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5S04_p, 1, function(x) {mean(x<=0.05)}))


LT4S1 <- read.csv("LGCP\\LGCP_1_TILING_4_SMALL.csv")[,-1]
LT4S08 <- read.csv("LGCP\\LGCP_08_TILING_4_SMALL.csv")[,-1]
LT4S06 <- read.csv("LGCP\\LGCP_06_TILING_4_SMALL.csv")[,-1]
LT4S04 <- read.csv("LGCP\\LGCP_04_TILING_4_SMALL.csv")[,-1]

LT4S1_p <- rbind(apply(LT4S1[2:1002,], 2, pval),
                apply(LT4S1[2006:3006,], 2, pval),
                apply(LT4S1[4010:5010,], 2, pval))
LT4S08_p <- rbind(apply(LT4S08[2:1002,], 2, pval),
                 apply(LT4S08[2006:3006,], 2, pval),
                 apply(LT4S08[4010:5010,], 2, pval))
LT4S06_p <- rbind(apply(LT4S06[2:1002,], 2, pval),
                 apply(LT4S06[2006:3006,], 2, pval),
                 apply(LT4S06[4010:5010,], 2, pval))
LT4S04_p <- rbind(apply(LT4S04[2:1002,], 2, pval),
                 apply(LT4S04[2006:3006,], 2, pval),
                 apply(LT4S04[4010:5010,], 2, pval))
POW_LT4S <- cbind(apply(LT4S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4S04_p, 1, function(x) {mean(x<=0.05)}))


LT3S1 <- read.csv("LGCP\\LGCP_1_TILING_3_SMALL.csv")[,-1]
LT3S08 <- read.csv("LGCP\\LGCP_08_TILING_3_SMALL.csv")[,-1]
LT3S06 <- read.csv("LGCP\\LGCP_06_TILING_3_SMALL.csv")[,-1]
LT3S04 <- read.csv("LGCP\\LGCP_04_TILING_3_SMALL.csv")[,-1]

LT3S1_p <- rbind(apply(LT3S1[2:1002,], 2, pval),
                apply(LT3S1[2006:3006,], 2, pval),
                apply(LT3S1[4010:5010,], 2, pval))
LT3S08_p <- rbind(apply(LT3S08[2:1002,], 2, pval),
                 apply(LT3S08[2006:3006,], 2, pval),
                 apply(LT3S08[4010:5010,], 2, pval))
LT3S06_p <- rbind(apply(LT3S06[2:1002,], 2, pval),
                 apply(LT3S06[2006:3006,], 2, pval),
                 apply(LT3S06[4010:5010,], 2, pval))
LT3S04_p <- rbind(apply(LT3S04[2:1002,], 2, pval),
                 apply(LT3S04[2006:3006,], 2, pval),
                 apply(LT3S04[4010:5010,], 2, pval))
POW_LT3S <- cbind(apply(LT3S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT3S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT3S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT3S04_p, 1, function(x) {mean(x<=0.05)}))


LT2S1 <- read.csv("LGCP\\LGCP_1_TILING_2_SMALL.csv")[,-1]
LT2S08 <- read.csv("LGCP\\LGCP_08_TILING_2_SMALL.csv")[,-1]
LT2S06 <- read.csv("LGCP\\LGCP_06_TILING_2_SMALL.csv")[,-1]
LT2S04 <- read.csv("LGCP\\LGCP_04_TILING_2_SMALL.csv")[,-1]

LT2S1_p <- rbind(apply(LT2S1[2:1002,], 2, pval),
                apply(LT2S1[2006:3006,], 2, pval),
                apply(LT2S1[4010:5010,], 2, pval))
LT2S08_p <- rbind(apply(LT2S08[2:1002,], 2, pval),
                 apply(LT2S08[2006:3006,], 2, pval),
                 apply(LT2S08[4010:5010,], 2, pval))
LT2S06_p <- rbind(apply(LT2S06[2:1002,], 2, pval),
                 apply(LT2S06[2006:3006,], 2, pval),
                 apply(LT2S06[4010:5010,], 2, pval))
LT2S04_p <- rbind(apply(LT2S04[2:1002,], 2, pval),
                 apply(LT2S04[2006:3006,], 2, pval),
                 apply(LT2S04[4010:5010,], 2, pval))
POW_LT2S <- cbind(apply(LT2S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT2S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT2S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT2S04_p, 1, function(x) {mean(x<=0.05)}))


LLB1 <- read.csv("LGCP\\LGCP_LGCP_1_BIG.csv")[,-1]
LLB08 <- read.csv("LGCP\\LGCP_LGCP_08_BIG.csv")[,-1]
LLB06 <- read.csv("LGCP\\LGCP_LGCP_06_BIG.csv")[,-1]
LLB04 <- read.csv("LGCP\\LGCP_LGCP_04_BIG.csv")[,-1]

LLB1_p <- rbind(apply(LLB1[2:1002,], 2, pval),
                apply(LLB1[2006:3006,], 2, pval),
                apply(LLB1[4010:5010,], 2, pval))
LLB08_p <- rbind(apply(LLB08[2:1002,], 2, pval),
                 apply(LLB08[2006:3006,], 2, pval),
                 apply(LLB08[4010:5010,], 2, pval))
LLB06_p <- rbind(apply(LLB06[2:1002,], 2, pval),
                 apply(LLB06[2006:3006,], 2, pval),
                 apply(LLB06[4010:5010,], 2, pval))
LLB04_p <- rbind(apply(LLB04[2:1002,], 2, pval),
                 apply(LLB04[2006:3006,], 2, pval),
                 apply(LLB04[4010:5010,], 2, pval))
POW_LLB <- cbind(apply(LLB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LLB04_p, 1, function(x) {mean(x<=0.05)}))


LTB1 <- read.csv("LGCP\\LGCP_THOMAS_1_BIG_long_short.csv")[,-1]
LTB08 <- read.csv("LGCP\\LGCP_THOMAS_08_BIG_long_short.csv")[,-1]
LTB06 <- read.csv("LGCP\\LGCP_THOMAS_06_BIG_long_short.csv")[,-1]
LTB04 <- read.csv("LGCP\\LGCP_THOMAS_04_BIG_long_short.csv")[,-1]

LTB1_p <- rbind(apply(LTB1[2:1002,], 2, pval),
                apply(LTB1[2006:3006,], 2, pval),
                apply(LTB1[4010:5010,], 2, pval))
LTB08_p <- rbind(apply(LTB08[2:1002,], 2, pval),
                 apply(LTB08[2006:3006,], 2, pval),
                 apply(LTB08[4010:5010,], 2, pval))
LTB06_p <- rbind(apply(LTB06[2:1002,], 2, pval),
                 apply(LTB06[2006:3006,], 2, pval),
                 apply(LTB06[4010:5010,], 2, pval))
LTB04_p <- rbind(apply(LTB04[2:1002,], 2, pval),
                 apply(LTB04[2006:3006,], 2, pval),
                 apply(LTB04[4010:5010,], 2, pval))
POW_LTB <- cbind(apply(LTB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LTB04_p, 1, function(x) {mean(x<=0.05)}))


LSB1 <- read.csv("LGCP\\LGCP_SR_1_BIG.csv")[,-1]
LSB08 <- read.csv("LGCP\\LGCP_SR_08_BIG.csv")[,-1]
LSB06 <- read.csv("LGCP\\LGCP_SR_06_BIG.csv")[,-1]
LSB04 <- read.csv("LGCP\\LGCP_SR_04_BIG.csv")[,-1]

LSB1_p <- rbind(apply(LSB1[2:101,], 2, pval_SR),
                apply(LSB1[204:305,], 2, pval_SR),
                apply(LSB1[406:505,], 2, pval_SR))
LSB08_p <- rbind(apply(LSB08[2:101,], 2, pval_SR),
                 apply(LSB08[204:305,], 2, pval_SR),
                 apply(LSB08[406:505,], 2, pval_SR))
LSB06_p <- rbind(apply(LSB06[2:101,], 2, pval_SR),
                 apply(LSB06[204:305,], 2, pval_SR),
                 apply(LSB06[406:505,], 2, pval_SR))
LSB04_p <- rbind(apply(LSB04[2:101,], 2, pval_SR),
                 apply(LSB04[204:305,], 2, pval_SR),
                 apply(LSB04[406:505,], 2, pval_SR))
POW_LSB <- cbind(apply(LSB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LSB04_p, 1, function(x) {mean(x<=0.05)}))


LT8B1 <- read.csv("LGCP\\LGCP_1_TILING_8_BIG.csv")[,-1]
LT8B08 <- read.csv("LGCP\\LGCP_08_TILING_8_BIG.csv")[,-1]
LT8B06 <- read.csv("LGCP\\LGCP_06_TILING_8_BIG.csv")[,-1]
LT8B04 <- read.csv("LGCP\\LGCP_04_TILING_8_BIG.csv")[,-1]

LT8B1_p <- rbind(apply(LT8B1[2:1002,], 2, pval),
                apply(LT8B1[2006:3006,], 2, pval),
                apply(LT8B1[4010:5010,], 2, pval))
LT8B08_p <- rbind(apply(LT8B08[2:1002,], 2, pval),
                 apply(LT8B08[2006:3006,], 2, pval),
                 apply(LT8B08[4010:5010,], 2, pval))
LT8B06_p <- rbind(apply(LT8B06[2:1002,], 2, pval),
                 apply(LT8B06[2006:3006,], 2, pval),
                 apply(LT8B06[4010:5010,], 2, pval))
LT8B04_p <- rbind(apply(LT8B04[2:1002,], 2, pval),
                 apply(LT8B04[2006:3006,], 2, pval),
                 apply(LT8B04[4010:5010,], 2, pval))
POW_LT8B <- cbind(apply(LT8B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT8B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT8B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT8B04_p, 1, function(x) {mean(x<=0.05)}))


LT6B1 <- read.csv("LGCP\\LGCP_1_TILING_6_BIG.csv")[,-1]
LT6B08 <- read.csv("LGCP\\LGCP_08_TILING_6_BIG.csv")[,-1]
LT6B06 <- read.csv("LGCP\\LGCP_06_TILING_6_BIG.csv")[,-1]
LT6B04 <- read.csv("LGCP\\LGCP_04_TILING_6_BIG.csv")[,-1]

LT6B1_p <- rbind(apply(LT6B1[2:1002,], 2, pval),
                apply(LT6B1[2006:3006,], 2, pval),
                apply(LT6B1[4010:5010,], 2, pval))
LT6B08_p <- rbind(apply(LT6B08[2:1002,], 2, pval),
                 apply(LT6B08[2006:3006,], 2, pval),
                 apply(LT6B08[4010:5010,], 2, pval))
LT6B06_p <- rbind(apply(LT6B06[2:1002,], 2, pval),
                 apply(LT6B06[2006:3006,], 2, pval),
                 apply(LT6B06[4010:5010,], 2, pval))
LT6B04_p <- rbind(apply(LT6B04[2:1002,], 2, pval),
                 apply(LT6B04[2006:3006,], 2, pval),
                 apply(LT6B04[4010:5010,], 2, pval))
POW_LT6B <- cbind(apply(LT6B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT6B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT6B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT6B04_p, 1, function(x) {mean(x<=0.05)}))


LT5B1 <- read.csv("LGCP\\LGCP_1_TILING_5_BIG.csv")[,-1]
LT5B08 <- read.csv("LGCP\\LGCP_08_TILING_5_BIG.csv")[,-1]
LT5B06 <- read.csv("LGCP\\LGCP_06_TILING_5_BIG.csv")[,-1]
LT5B04 <- read.csv("LGCP\\LGCP_04_TILING_5_BIG.csv")[,-1]

LT5B1_p <- rbind(apply(LT5B1[2:1002,], 2, pval),
                apply(LT5B1[2006:3006,], 2, pval),
                apply(LT5B1[4010:5010,], 2, pval))
LT5B08_p <- rbind(apply(LT5B08[2:1002,], 2, pval),
                 apply(LT5B08[2006:3006,], 2, pval),
                 apply(LT5B08[4010:5010,], 2, pval))
LT5B06_p <- rbind(apply(LT5B06[2:1002,], 2, pval),
                 apply(LT5B06[2006:3006,], 2, pval),
                 apply(LT5B06[4010:5010,], 2, pval))
LT5B04_p <- rbind(apply(LT5B04[2:1002,], 2, pval),
                 apply(LT5B04[2006:3006,], 2, pval),
                 apply(LT5B04[4010:5010,], 2, pval))
POW_LT5B <- cbind(apply(LT5B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT5B04_p, 1, function(x) {mean(x<=0.05)}))


LT4B1 <- read.csv("LGCP\\LGCP_1_TILING_4_BIG.csv")[,-1]
LT4B08 <- read.csv("LGCP\\LGCP_08_TILING_4_BIG.csv")[,-1]
LT4B06 <- read.csv("LGCP\\LGCP_06_TILING_4_BIG.csv")[,-1]
LT4B04 <- read.csv("LGCP\\LGCP_04_TILING_4_BIG.csv")[,-1]

LT4B1_p <- rbind(apply(LT4B1[2:1002,], 2, pval),
                apply(LT4B1[2006:3006,], 2, pval),
                apply(LT4B1[4010:5010,], 2, pval))
LT4B08_p <- rbind(apply(LT4B08[2:1002,], 2, pval),
                 apply(LT4B08[2006:3006,], 2, pval),
                 apply(LT4B08[4010:5010,], 2, pval))
LT4B06_p <- rbind(apply(LT4B06[2:1002,], 2, pval),
                 apply(LT4B06[2006:3006,], 2, pval),
                 apply(LT4B06[4010:5010,], 2, pval))
LT4B04_p <- rbind(apply(LT4B04[2:1002,], 2, pval),
                 apply(LT4B04[2006:3006,], 2, pval),
                 apply(LT4B04[4010:5010,], 2, pval))
POW_LT4B <- cbind(apply(LT4B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LT4B04_p, 1, function(x) {mean(x<=0.05)}))


LOB1 <- read.csv("LGCP\\LGCP_ORACLE_1_BIG.csv")[,-1]
LOB08 <- read.csv("LGCP\\LGCP_ORACLE_08_BIG.csv")[,-1]
LOB06 <- read.csv("LGCP\\LGCP_ORACLE_06_BIG.csv")[,-1]
LOB04 <- read.csv("LGCP\\LGCP_ORACLE_04_BIG.csv")[,-1]

LOB1_p <- rbind(apply(LOB1[2:1002,], 2, pval),
                apply(LOB1[2006:3006,], 2, pval),
                apply(LOB1[4010:5010,], 2, pval))
LOB08_p <- rbind(apply(LOB08[2:1002,], 2, pval),
                 apply(LOB08[2006:3006,], 2, pval),
                 apply(LOB08[4010:5010,], 2, pval))
LOB06_p <- rbind(apply(LOB06[2:1002,], 2, pval),
                 apply(LOB06[2006:3006,], 2, pval),
                 apply(LOB06[4010:5010,], 2, pval))
LOB04_p <- rbind(apply(LOB04[2:1002,], 2, pval),
                 apply(LOB04[2006:3006,], 2, pval),
                 apply(LOB04[4010:5010,], 2, pval))
POW_LOB <- cbind(apply(LOB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOB04_p, 1, function(x) {mean(x<=0.05)}))

LOS1 <- read.csv("LGCP\\LGCP_ORACLE_1_SMALL.csv")[,-1]
LOS08 <- read.csv("LGCP\\LGCP_ORACLE_08_SMALL.csv")[,-1]
LOS06 <- read.csv("LGCP\\LGCP_ORACLE_06_SMALL.csv")[,-1]
LOS04 <- read.csv("LGCP\\LGCP_ORACLE_04_SMALL.csv")[,-1]

LOS1_p <- rbind(apply(LOS1[2:1002,], 2, pval),
                apply(LOS1[2006:3006,], 2, pval),
                apply(LOS1[4010:5010,], 2, pval))
LOS08_p <- rbind(apply(LOS08[2:1002,], 2, pval),
                 apply(LOS08[2006:3006,], 2, pval),
                 apply(LOS08[4010:5010,], 2, pval))
LOS06_p <- rbind(apply(LOS06[2:1002,], 2, pval),
                 apply(LOS06[2006:3006,], 2, pval),
                 apply(LOS06[4010:5010,], 2, pval))
LOS04_p <- rbind(apply(LOS04[2:1002,], 2, pval),
                 apply(LOS04[2006:3006,], 2, pval),
                 apply(LOS04[4010:5010,], 2, pval))
POW_LOS <- cbind(apply(LOS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(LOS04_p, 1, function(x) {mean(x<=0.05)}))


POW_LS <- rbind(POW_LOS, POW_LLS,POW_LTS, POW_LT2S, POW_LT3S, POW_LT4S, POW_LT5S, POW_LSS)
colnames(POW_LS) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_LS) <- c("LOS_G", "LOS_K", "LOS_T",
                      "LLS_G", "LLS_K", "LLS_T",
                      "LTS_G", "LTS_K", "LTS_T",
                      "LT2S_G", "LT2S_K", "LT2S_T",
                      "LT3S_G", "LT3S_K", "LT3S_T",
                      "LT4S_G", "LT4S_K", "LT4S_T",
                      "LT5S_G", "LT5S_K", "LT5S_T",
                      "LSS_G", "LSS_K", "LSS_T")

write.csv(POW_LS, "POW_LS.csv")


POW_LB <- rbind(POW_LOB, POW_LLB,POW_LTB, POW_LT4B, POW_LT5B, POW_LT6B, POW_LT8B, POW_LSB)
colnames(POW_LB) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_LB) <- c("LOB_G", "LOB_K", "LOB_T",
                      "LLB_G", "LLB_K", "LLB_T",
                      "LTB_G", "LTB_K", "LTB_T",
                      "LT4B_G", "LT4B_K", "LT4B_T",
                      "LT5B_G", "LT5B_K", "LT5B_T",
                      "LT6B_G", "LT6B_K", "LT6B_T",
                      "LT8B_G", "LT8B_K", "LT8B_T",
                      "LSB_G", "LSB_K", "LSB_T")

write.csv(POW_LB, "POW_LB.csv")


##### 
#GIBBS


GSS1 <- read.csv("GIBBS\\GIBBS_STRAUSS_1_SMALL.csv")[,-1]
GSS08 <- read.csv("GIBBS\\GIBBS_STRAUSS_08_SMALL.csv")[,-1]
GSS06 <- read.csv("GIBBS\\GIBBS_STRAUSS_06_SMALL.csv")[,-1]
GSS04 <- read.csv("GIBBS\\GIBBS_STRAUSS_04_SMALL.csv")[,-1]

GSS1_p <- rbind(apply(GSS1[2:1002,], 2, pval),
                apply(GSS1[2006:3006,], 2, pval),
                apply(GSS1[4010:5010,], 2, pval))
GSS08_p <- rbind(apply(GSS08[2:1002,], 2, pval),
                 apply(GSS08[2006:3006,], 2, pval),
                 apply(GSS08[4010:5010,], 2, pval))
GSS06_p <- rbind(apply(GSS06[2:1002,], 2, pval),
                 apply(GSS06[2006:3006,], 2, pval),
                 apply(GSS06[4010:5010,], 2, pval))
GSS04_p <- rbind(apply(GSS04[2:1002,], 2, pval),
                 apply(GSS04[2006:3006,], 2, pval),
                 apply(GSS04[4010:5010,], 2, pval))
POW_GSS <- cbind(apply(GSS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSS04_p, 1, function(x) {mean(x<=0.05)}))



GSRS1 <- read.csv("GIBBS\\GIBBS_SR_1_SMALL.csv")[,-1]
GSRS08 <- read.csv("GIBBS\\GIBBS_SR_08_SMALL.csv")[,-1]
GSRS06 <- read.csv("GIBBS\\GIBBS_SR_06_SMALL.csv")[,-1]
GSRS04 <- read.csv("GIBBS\\GIBBS_SR_04_SMALL.csv")[,-1]

GSRS1_p <- rbind(apply(GSRS1[2:101,], 2, pval_SR),
                apply(GSRS1[204:305,], 2, pval_SR),
                apply(GSRS1[406:505,], 2, pval_SR))
GSRS08_p <- rbind(apply(GSRS08[2:101,], 2, pval_SR),
                 apply(GSRS08[204:305,], 2, pval_SR),
                 apply(GSRS08[406:505,], 2, pval_SR))
GSRS06_p <- rbind(apply(GSRS06[2:101,], 2, pval_SR),
                 apply(GSRS06[204:305,], 2, pval_SR),
                 apply(GSRS06[406:505,], 2, pval_SR))
GSRS04_p <- rbind(apply(GSRS04[2:101,], 2, pval_SR),
                 apply(GSRS04[204:305,], 2, pval_SR),
                 apply(GSRS04[406:505,], 2, pval_SR))
POW_GSRS <- cbind(apply(GSRS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSRS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSRS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSRS04_p, 1, function(x) {mean(x<=0.05)}))

GT5S1 <- read.csv("GIBBS\\GIBBS_1_TILING_5_SMALL.csv")[,-1]
GT5S08 <- read.csv("GIBBS\\GIBBS_08_TILING_5_SMALL.csv")[,-1]
GT5S06 <- read.csv("GIBBS\\GIBBS_06_TILING_5_SMALL.csv")[,-1]
GT5S04 <- read.csv("GIBBS\\GIBBS_04_TILING_5_SMALL.csv")[,-1]

GT5S1_p <- rbind(apply(GT5S1[2:1002,], 2, pval),
                apply(GT5S1[2006:3006,], 2, pval),
                apply(GT5S1[4010:5010,], 2, pval))
GT5S08_p <- rbind(apply(GT5S08[2:1002,], 2, pval),
                 apply(GT5S08[2006:3006,], 2, pval),
                 apply(GT5S08[4010:5010,], 2, pval))
GT5S06_p <- rbind(apply(GT5S06[2:1002,], 2, pval),
                 apply(GT5S06[2006:3006,], 2, pval),
                 apply(GT5S06[4010:5010,], 2, pval))
GT5S04_p <- rbind(apply(GT5S04[2:1002,], 2, pval),
                 apply(GT5S04[2006:3006,], 2, pval),
                 apply(GT5S04[4010:5010,], 2, pval))
POW_GT5S <- cbind(apply(GT5S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5S04_p, 1, function(x) {mean(x<=0.05)}))

GT4S1 <- read.csv("GIBBS\\GIBBS_1_TILING_4_SMALL.csv")[,-1]
GT4S08 <- read.csv("GIBBS\\GIBBS_08_TILING_4_SMALL.csv")[,-1]
GT4S06 <- read.csv("GIBBS\\GIBBS_06_TILING_4_SMALL.csv")[,-1]
GT4S04 <- read.csv("GIBBS\\GIBBS_04_TILING_4_SMALL.csv")[,-1]

GT4S1_p <- rbind(apply(GT4S1[2:1002,], 2, pval),
                apply(GT4S1[2006:3006,], 2, pval),
                apply(GT4S1[4010:5010,], 2, pval))
GT4S08_p <- rbind(apply(GT4S08[2:1002,], 2, pval),
                 apply(GT4S08[2006:3006,], 2, pval),
                 apply(GT4S08[4010:5010,], 2, pval))
GT4S06_p <- rbind(apply(GT4S06[2:1002,], 2, pval),
                 apply(GT4S06[2006:3006,], 2, pval),
                 apply(GT4S06[4010:5010,], 2, pval))
GT4S04_p <- rbind(apply(GT4S04[2:1002,], 2, pval),
                 apply(GT4S04[2006:3006,], 2, pval),
                 apply(GT4S04[4010:5010,], 2, pval))
POW_GT4S <- cbind(apply(GT4S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4S04_p, 1, function(x) {mean(x<=0.05)}))

GT3S1 <- read.csv("GIBBS\\GIBBS_1_TILING_3_SMALL.csv")[,-1]
GT3S08 <- read.csv("GIBBS\\GIBBS_08_TILING_3_SMALL.csv")[,-1]
GT3S06 <- read.csv("GIBBS\\GIBBS_06_TILING_3_SMALL.csv")[,-1]
GT3S04 <- read.csv("GIBBS\\GIBBS_04_TILING_3_SMALL.csv")[,-1]

GT3S1_p <- rbind(apply(GT3S1[2:1002,], 2, pval),
                apply(GT3S1[2006:3006,], 2, pval),
                apply(GT3S1[4010:5010,], 2, pval))
GT3S08_p <- rbind(apply(GT3S08[2:1002,], 2, pval),
                 apply(GT3S08[2006:3006,], 2, pval),
                 apply(GT3S08[4010:5010,], 2, pval))
GT3S06_p <- rbind(apply(GT3S06[2:1002,], 2, pval),
                 apply(GT3S06[2006:3006,], 2, pval),
                 apply(GT3S06[4010:5010,], 2, pval))
GT3S04_p <- rbind(apply(GT3S04[2:1002,], 2, pval),
                 apply(GT3S04[2006:3006,], 2, pval),
                 apply(GT3S04[4010:5010,], 2, pval))
POW_GT3S <- cbind(apply(GT3S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT3S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT3S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT3S04_p, 1, function(x) {mean(x<=0.05)}))

GT2S1 <- read.csv("GIBBS\\GIBBS_1_TILING_2_SMALL.csv")[,-1]
GT2S08 <- read.csv("GIBBS\\GIBBS_08_TILING_2_SMALL.csv")[,-1]
GT2S06 <- read.csv("GIBBS\\GIBBS_06_TILING_2_SMALL.csv")[,-1]
GT2S04 <- read.csv("GIBBS\\GIBBS_04_TILING_2_SMALL.csv")[,-1]

GT2S1_p <- rbind(apply(GT2S1[2:1002,], 2, pval),
                apply(GT2S1[2006:3006,], 2, pval),
                apply(GT2S1[4010:5010,], 2, pval))
GT2S08_p <- rbind(apply(GT2S08[2:1002,], 2, pval),
                 apply(GT2S08[2006:3006,], 2, pval),
                 apply(GT2S08[4010:5010,], 2, pval))
GT2S06_p <- rbind(apply(GT2S06[2:1002,], 2, pval),
                 apply(GT2S06[2006:3006,], 2, pval),
                 apply(GT2S06[4010:5010,], 2, pval))
GT2S04_p <- rbind(apply(GT2S04[2:1002,], 2, pval),
                 apply(GT2S04[2006:3006,], 2, pval),
                 apply(GT2S04[4010:5010,], 2, pval))
POW_GT2S <- cbind(apply(GT2S1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT2S08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT2S06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT2S04_p, 1, function(x) {mean(x<=0.05)}))

GOS1 <- read.csv("GIBBS\\GIBBS_ORACLE_1_SMALL.csv")[,-1]
GOS08 <- read.csv("GIBBS\\GIBBS_ORACLE_08_SMALL.csv")[,-1]
GOS06 <- read.csv("GIBBS\\GIBBS_ORACLE_06_SMALL.csv")[,-1]
GOS04 <- read.csv("GIBBS\\GIBBS_ORACLE_04_SMALL.csv")[,-1]

GOS1_p <- rbind(apply(GOS1[2:1002,], 2, pval),
                apply(GOS1[2006:3006,], 2, pval),
                apply(GOS1[4010:5010,], 2, pval))
GOS08_p <- rbind(apply(GOS08[2:1002,], 2, pval),
                 apply(GOS08[2006:3006,], 2, pval),
                 apply(GOS08[4010:5010,], 2, pval))
GOS06_p <- rbind(apply(GOS06[2:1002,], 2, pval),
                 apply(GOS06[2006:3006,], 2, pval),
                 apply(GOS06[4010:5010,], 2, pval))
GOS04_p <- rbind(apply(GOS04[2:1002,], 2, pval),
                 apply(GOS04[2006:3006,], 2, pval),
                 apply(GOS04[4010:5010,], 2, pval))
POW_GOS <- cbind(apply(GOS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOS04_p, 1, function(x) {mean(x<=0.05)}))

GSB1 <- read.csv("GIBBS\\GIBBS_STRAUSS_1_BIG.csv")[,-1]
GSB08 <- read.csv("GIBBS\\GIBBS_STRAUSS_08_BIG.csv")[,-1]
GSB06 <- read.csv("GIBBS\\GIBBS_STRAUSS_06_BIG.csv")[,-1]
GSB04 <- read.csv("GIBBS\\GIBBS_STRAUSS_04_BIG.csv")[,-1]

GSB1_p <- rbind(apply(GSB1[2:1002,], 2, pval),
                apply(GSB1[2006:3006,], 2, pval),
                apply(GSB1[4010:5010,], 2, pval))
GSB08_p <- rbind(apply(GSB08[2:1002,], 2, pval),
                 apply(GSB08[2006:3006,], 2, pval),
                 apply(GSB08[4010:5010,], 2, pval))
GSB06_p <- rbind(apply(GSB06[2:1002,], 2, pval),
                 apply(GSB06[2006:3006,], 2, pval),
                 apply(GSB06[4010:5010,], 2, pval))
GSB04_p <- rbind(apply(GSB04[2:1002,], 2, pval),
                 apply(GSB04[2006:3006,], 2, pval),
                 apply(GSB04[4010:5010,], 2, pval))
POW_GSB <- cbind(apply(GSB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GSB04_p, 1, function(x) {mean(x<=0.05)}))

GSRB1 <- read.csv("GIBBS\\GIBBS_SR_1_BIG.csv")[,-1]
GSRB08 <- read.csv("GIBBS\\GIBBS_SR_08_BIG.csv")[,-1]
GSRB06 <- read.csv("GIBBS\\GIBBS_SR_06_BIG.csv")[,-1]
GSRB04 <- read.csv("GIBBS\\GIBBS_SR_04_BIG.csv")[,-1]

GSRB1_p <- rbind(apply(GSRB1[2:101,], 2, pval_SR),
                 apply(GSRB1[204:305,], 2, pval_SR),
                 apply(GSRB1[406:505,], 2, pval_SR))
GSRB08_p <- rbind(apply(GSRB08[2:101,], 2, pval_SR),
                  apply(GSRB08[204:305,], 2, pval_SR),
                  apply(GSRB08[406:505,], 2, pval_SR))
GSRB06_p <- rbind(apply(GSRB06[2:101,], 2, pval_SR),
                  apply(GSRB06[204:305,], 2, pval_SR),
                  apply(GSRB06[406:505,], 2, pval_SR))
GSRB04_p <- rbind(apply(GSRB04[2:101,], 2, pval_SR),
                  apply(GSRB04[204:305,], 2, pval_SR),
                  apply(GSRB04[406:505,], 2, pval_SR))
POW_GSRB <- cbind(apply(GSRB1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(GSRB08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(GSRB06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(GSRB04_p, 1, function(x) {mean(x<=0.05)}))

GT8B1 <- read.csv("GIBBS\\GIBBS_1_TILING_8_BIG.csv")[,-1]
GT8B08 <- read.csv("GIBBS\\GIBBS_08_TILING_8_BIG.csv")[,-1]
GT8B06 <- read.csv("GIBBS\\GIBBS_06_TILING_8_BIG.csv")[,-1]
GT8B04 <- read.csv("GIBBS\\GIBBS_04_TILING_8_BIG.csv")[,-1]

GT8B1_p <- rbind(apply(GT8B1[2:1002,], 2, pval),
                apply(GT8B1[2006:3006,], 2, pval),
                apply(GT8B1[4010:5010,], 2, pval))
GT8B08_p <- rbind(apply(GT8B08[2:1002,], 2, pval),
                 apply(GT8B08[2006:3006,], 2, pval),
                 apply(GT8B08[4010:5010,], 2, pval))
GT8B06_p <- rbind(apply(GT8B06[2:1002,], 2, pval),
                 apply(GT8B06[2006:3006,], 2, pval),
                 apply(GT8B06[4010:5010,], 2, pval))
GT8B04_p <- rbind(apply(GT8B04[2:1002,], 2, pval),
                 apply(GT8B04[2006:3006,], 2, pval),
                 apply(GT8B04[4010:5010,], 2, pval))
POW_GT8B <- cbind(apply(GT8B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT8B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT8B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT8B04_p, 1, function(x) {mean(x<=0.05)}))

GT6B1 <- read.csv("GIBBS\\GIBBS_1_TILING_6_BIG.csv")[,-1]
GT6B08 <- read.csv("GIBBS\\GIBBS_08_TILING_6_BIG.csv")[,-1]
GT6B06 <- read.csv("GIBBS\\GIBBS_06_TILING_6_BIG.csv")[,-1]
GT6B04 <- read.csv("GIBBS\\GIBBS_04_TILING_6_BIG.csv")[,-1]

GT6B1_p <- rbind(apply(GT6B1[2:1002,], 2, pval),
                apply(GT6B1[2006:3006,], 2, pval),
                apply(GT6B1[4010:5010,], 2, pval))
GT6B08_p <- rbind(apply(GT6B08[2:1002,], 2, pval),
                 apply(GT6B08[2006:3006,], 2, pval),
                 apply(GT6B08[4010:5010,], 2, pval))
GT6B06_p <- rbind(apply(GT6B06[2:1002,], 2, pval),
                 apply(GT6B06[2006:3006,], 2, pval),
                 apply(GT6B06[4010:5010,], 2, pval))
GT6B04_p <- rbind(apply(GT6B04[2:1002,], 2, pval),
                 apply(GT6B04[2006:3006,], 2, pval),
                 apply(GT6B04[4010:5010,], 2, pval))
POW_GT6B <- cbind(apply(GT6B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT6B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT6B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT6B04_p, 1, function(x) {mean(x<=0.05)}))

GT4B1 <- read.csv("GIBBS\\GIBBS_1_TILING_4_BIG.csv")[,-1]
GT4B08 <- read.csv("GIBBS\\GIBBS_08_TILING_4_BIG.csv")[,-1]
GT4B06 <- read.csv("GIBBS\\GIBBS_06_TILING_4_BIG.csv")[,-1]
GT4B04 <- read.csv("GIBBS\\GIBBS_04_TILING_4_BIG.csv")[,-1]

GT4B1_p <- rbind(apply(GT4B1[2:1002,], 2, pval),
                apply(GT4B1[2006:3006,], 2, pval),
                apply(GT4B1[4010:5010,], 2, pval))
GT4B08_p <- rbind(apply(GT4B08[2:1002,], 2, pval),
                 apply(GT4B08[2006:3006,], 2, pval),
                 apply(GT4B08[4010:5010,], 2, pval))
GT4B06_p <- rbind(apply(GT4B06[2:1002,], 2, pval),
                 apply(GT4B06[2006:3006,], 2, pval),
                 apply(GT4B06[4010:5010,], 2, pval))
GT4B04_p <- rbind(apply(GT4B04[2:1002,], 2, pval),
                 apply(GT4B04[2006:3006,], 2, pval),
                 apply(GT4B04[4010:5010,], 2, pval))
POW_GT4B <- cbind(apply(GT4B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT4B04_p, 1, function(x) {mean(x<=0.05)}))

GT5B1 <- read.csv("GIBBS\\GIBBS_1_TILING_5_BIG.csv")[,-1]
GT5B08 <- read.csv("GIBBS\\GIBBS_08_TILING_5_BIG.csv")[,-1]
GT5B06 <- read.csv("GIBBS\\GIBBS_06_TILING_5_BIG.csv")[,-1]
GT5B04 <- read.csv("GIBBS\\GIBBS_04_TILING_5_BIG.csv")[,-1]

GT5B1_p <- rbind(apply(GT5B1[2:1002,], 2, pval),
                apply(GT5B1[2006:3006,], 2, pval),
                apply(GT5B1[4010:5010,], 2, pval))
GT5B08_p <- rbind(apply(GT5B08[2:1002,], 2, pval),
                 apply(GT5B08[2006:3006,], 2, pval),
                 apply(GT5B08[4010:5010,], 2, pval))
GT5B06_p <- rbind(apply(GT5B06[2:1002,], 2, pval),
                 apply(GT5B06[2006:3006,], 2, pval),
                 apply(GT5B06[4010:5010,], 2, pval))
GT5B04_p <- rbind(apply(GT5B04[2:1002,], 2, pval),
                 apply(GT5B04[2006:3006,], 2, pval),
                 apply(GT5B04[4010:5010,], 2, pval))
POW_GT5B <- cbind(apply(GT5B1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5B08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5B06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GT5B04_p, 1, function(x) {mean(x<=0.05)}))

GOB1 <- read.csv("GIBBS\\GIBBS_ORACLE_1_BIG.csv")[,-1]
GOB08 <- read.csv("GIBBS\\GIBBS_ORACLE_08_BIG.csv")[,-1]
GOB06 <- read.csv("GIBBS\\GIBBS_ORACLE_06_BIG.csv")[,-1]
GOB04 <- read.csv("GIBBS\\GIBBS_ORACLE_04_BIG.csv")[,-1]

GOB1_p <- rbind(apply(GOB1[2:1002,], 2, pval),
                apply(GOB1[2006:3006,], 2, pval),
                apply(GOB1[4010:5010,], 2, pval))
GOB08_p <- rbind(apply(GOB08[2:1002,], 2, pval),
                 apply(GOB08[2006:3006,], 2, pval),
                 apply(GOB08[4010:5010,], 2, pval))
GOB06_p <- rbind(apply(GOB06[2:1002,], 2, pval),
                 apply(GOB06[2006:3006,], 2, pval),
                 apply(GOB06[4010:5010,], 2, pval))
GOB04_p <- rbind(apply(GOB04[2:1002,], 2, pval),
                 apply(GOB04[2006:3006,], 2, pval),
                 apply(GOB04[4010:5010,], 2, pval))
POW_GOB <- cbind(apply(GOB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(GOB04_p, 1, function(x) {mean(x<=0.05)}))

POW_GS <- rbind(POW_GOS, POW_GSS, POW_GT2S, POW_GT3S, POW_GT4S, POW_GT5S, POW_GSRS)
colnames(POW_GS) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_GS) <- c("GOS_G", "GOS_K", "GOS_T",
                      "GSS_G", "GSS_K", "GSS_T",
                      "GT2S_G", "GT2S_K", "GT2S_T",
                      "GT3S_G", "GT3S_K", "GT3S_T",
                      "GT4S_G", "GT4S_K", "GT4S_T",
                      "GT5S_G", "GT5S_K", "GT5S_T",
                      "GSRS_G", "GSRS_K", "GSRS_T")


write.csv(POW_GS, "POW_GS.csv")


POW_GB <- rbind(POW_GOB, POW_GSB, POW_GT4B, POW_GT5B, POW_GT6B, POW_GT8B, POW_GSRB)
colnames(POW_GB) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_GB) <- c("GOB_G", "GOB_K", "GOB_T",
                      "GSB_G", "GSB_K", "GSB_T",
                      "GT4B_G", "GT4B_K", "GT4B_T",
                      "GT5B_G", "GT5B_K", "GT5B_T",
                      "GT6B_G", "GT6B_K", "GT6B_T",
                      "GT8B_G", "GT8B_K", "GT8B_T",
                      "GSRB_G", "GSRB_K", "GSRB_T")
write.csv(POW_GB, "POW_GB.csv")



#########
# PLCP
PPS1 <- read.csv("PLCP\\PLCP_PLCP_1_SMALL.csv")[,-1]
PPS08 <- read.csv("PLCP\\PLCP_PLCP_08_SMALL.csv")[,-1]
PPS06 <- read.csv("PLCP\\PLCP_PLCP_06_SMALL.csv")[,-1]
PPS04 <- read.csv("PLCP\\PLCP_PLCP_04_SMALL.csv")[,-1]

PTS1 <- read.csv("PLCP\\PLCP_THOMAS_1_SMALL_long_short.csv")[,-1]
PTS08 <- read.csv("PLCP\\PLCP_THOMAS_08_SMALL_long_short.csv")[,-1]
PTS06 <- read.csv("PLCP\\PLCP_THOMAS_06_SMALL_long_short.csv")[,-1]
PTS04 <- read.csv("PLCP\\PLCP_THOMAS_04_SMALL_long_short.csv")[,-1]

PSS1 <- read.csv("PLCP\\PLCP_SR_1_SMALL.csv")[,-1]
PSS08 <- read.csv("PLCP\\PLCP_SR_08_SMALL.csv")[,-1]
PSS06 <- read.csv("PLCP\\PLCP_SR_06_SMALL.csv")[,-1]
PSS04 <- read.csv("PLCP\\PLCP_SR_04_SMALL.csv")[,-1]

PT5S1 <- read.csv("PLCP\\PLCP_1_TILING_5_SMALL.csv")[,-1]
PT5S08 <- read.csv("PLCP\\PLCP_08_TILING_5_SMALL.csv")[,-1]
PT5S06 <- read.csv("PLCP\\PLCP_06_TILING_5_SMALL.csv")[,-1]
PT5S04 <- read.csv("PLCP\\PLCP_04_TILING_5_SMALL.csv")[,-1]

PT4S1 <- read.csv("PLCP\\PLCP_1_TILING_4_SMALL.csv")[,-1]
PT4S08 <- read.csv("PLCP\\PLCP_08_TILING_4_SMALL.csv")[,-1]
PT4S06 <- read.csv("PLCP\\PLCP_06_TILING_4_SMALL.csv")[,-1]
PT4S04 <- read.csv("PLCP\\PLCP_04_TILING_4_SMALL.csv")[,-1]

PT3S1 <- read.csv("PLCP\\PLCP_1_TILING_3_SMALL.csv")[,-1]
PT3S08 <- read.csv("PLCP\\PLCP_08_TILING_3_SMALL.csv")[,-1]
PT3S06 <- read.csv("PLCP\\PLCP_06_TILING_3_SMALL.csv")[,-1]
PT3S04 <- read.csv("PLCP\\PLCP_04_TILING_3_SMALL.csv")[,-1]

PT2S1 <- read.csv("PLCP\\PLCP_1_TILING_2_SMALL.csv")[,-1]
PT2S08 <- read.csv("PLCP\\PLCP_08_TILING_2_SMALL.csv")[,-1]
PT2S06 <- read.csv("PLCP\\PLCP_06_TILING_2_SMALL.csv")[,-1]
PT2S04 <- read.csv("PLCP\\PLCP_04_TILING_2_SMALL.csv")[,-1]

POS1 <- read.csv("PLCP\\PLCP_ORACLE_1_SMALL.csv")[,-1]
POS08 <- read.csv("PLCP\\PLCP_ORACLE_08_SMALL.csv")[,-1]
POS06 <- read.csv("PLCP\\PLCP_ORACLE_06_SMALL.csv")[,-1]
POS04 <- read.csv("PLCP\\PLCP_ORACLE_04_SMALL.csv")[,-1]

PPB1 <- read.csv("PLCP\\PLCP_PLCP_1_BIG.csv")[,-1]
PPB08 <- read.csv("PLCP\\PLCP_PLCP_08_BIG.csv")[,-1]
PPB06 <- read.csv("PLCP\\PLCP_PLCP_06_BIG.csv")[,-1]
PPB04 <- read.csv("PLCP\\PLCP_PLCP_04_BIG.csv")[,-1]

PTB1 <- read.csv("PLCP\\PLCP_THOMAS_1_BIG_long_short.csv")[,-1]
PTB08 <- read.csv("PLCP\\PLCP_THOMAS_08_BIG_long_short.csv")[,-1]
PTB06 <- read.csv("PLCP\\PLCP_THOMAS_06_BIG_long_short.csv")[,-1]
PTB04 <- read.csv("PLCP\\PLCP_THOMAS_04_BIG_long_short.csv")[,-1]

PSB1 <- read.csv("PLCP\\PLCP_SR_1_BIG.csv")[,-1]
PSB08 <- read.csv("PLCP\\PLCP_SR_08_BIG.csv")[,-1]
PSB06 <- read.csv("PLCP\\PLCP_SR_06_BIG.csv")[,-1]
PSB04 <- read.csv("PLCP\\PLCP_SR_04_BIG.csv")[,-1]

PT8B1 <- read.csv("PLCP\\PLCP_1_TILING_8_BIG.csv")[,-1]
PT8B08 <- read.csv("PLCP\\PLCP_08_TILING_8_BIG.csv")[,-1]
PT8B06 <- read.csv("PLCP\\PLCP_06_TILING_8_BIG.csv")[,-1]
PT8B04 <- read.csv("PLCP\\PLCP_04_TILING_8_BIG.csv")[,-1]

PT6B1 <- read.csv("PLCP\\PLCP_1_TILING_6_BIG.csv")[,-1]
PT6B08 <- read.csv("PLCP\\PLCP_08_TILING_6_BIG.csv")[,-1]
PT6B06 <- read.csv("PLCP\\PLCP_06_TILING_6_BIG.csv")[,-1]
PT6B04 <- read.csv("PLCP\\PLCP_04_TILING_6_BIG.csv")[,-1]

PT5B1 <- read.csv("PLCP\\PLCP_1_TILING_5_BIG.csv")[,-1]
PT5B08 <- read.csv("PLCP\\PLCP_08_TILING_5_BIG.csv")[,-1]
PT5B06 <- read.csv("PLCP\\PLCP_06_TILING_5_BIG.csv")[,-1]
PT5B04 <- read.csv("PLCP\\PLCP_04_TILING_5_BIG.csv")[,-1]

PT4B1 <- read.csv("PLCP\\PLCP_1_TILING_4_BIG.csv")[,-1]
PT4B08 <- read.csv("PLCP\\PLCP_08_TILING_4_BIG.csv")[,-1]
PT4B06 <- read.csv("PLCP\\PLCP_06_TILING_4_BIG.csv")[,-1]
PT4B04 <- read.csv("PLCP\\PLCP_04_TILING_4_BIG.csv")[,-1]

POB1 <- read.csv("PLCP\\PLCP_ORACLE_1_BIG.csv")[,-1]
POB08 <- read.csv("PLCP\\PLCP_ORACLE_08_BIG.csv")[,-1]
POB06 <- read.csv("PLCP\\PLCP_ORACLE_06_BIG.csv")[,-1]
POB04 <- read.csv("PLCP\\PLCP_ORACLE_04_BIG.csv")[,-1]




PPS1_p <- rbind(apply(PPS1[2:1002,], 2, pval),
                apply(PPS1[2006:3006,], 2, pval),
                apply(PPS1[4010:5010,], 2, pval))
PPS08_p <- rbind(apply(PPS08[2:1002,], 2, pval),
                 apply(PPS08[2006:3006,], 2, pval),
                 apply(PPS08[4010:5010,], 2, pval))
PPS06_p <- rbind(apply(PPS06[2:1002,], 2, pval),
                 apply(PPS06[2006:3006,], 2, pval),
                 apply(PPS06[4010:5010,], 2, pval))
PPS04_p <- rbind(apply(PPS04[2:1002,], 2, pval),
                 apply(PPS04[2006:3006,], 2, pval),
                 apply(PPS04[4010:5010,], 2, pval))
POW_PLS <- cbind(apply(PPS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPS04_p, 1, function(x) {mean(x<=0.05)}))



PTS1_p <- rbind(apply(PTS1[2:1002,], 2, pval),
                apply(PTS1[2006:3006,], 2, pval),
                apply(PTS1[4010:5010,], 2, pval))
PTS08_p <- rbind(apply(PTS08[2:1002,], 2, pval),
                 apply(PTS08[2006:3006,], 2, pval),
                 apply(PTS08[4010:5010,], 2, pval))
PTS06_p <- rbind(apply(PTS06[2:1002,], 2, pval),
                 apply(PTS06[2006:3006,], 2, pval),
                 apply(PTS06[4010:5010,], 2, pval))
PTS04_p <- rbind(apply(PTS04[2:1002,], 2, pval),
                 apply(PTS04[2006:3006,], 2, pval),
                 apply(PTS04[4010:5010,], 2, pval))
POW_PTS <- cbind(apply(PTS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTS04_p, 1, function(x) {mean(x<=0.05)}))


PSS1_p <- rbind(apply(PSS1[2:101,], 2, pval_SR),
                apply(PSS1[204:305,], 2, pval_SR),
                apply(PSS1[406:505,], 2, pval_SR))
PSS08_p <- rbind(apply(PSS08[2:101,], 2, pval_SR),
                 apply(PSS08[204:305,], 2, pval_SR),
                 apply(PSS08[406:505,], 2, pval_SR))
PSS06_p <- rbind(apply(PSS06[2:101,], 2, pval_SR),
                 apply(PSS06[204:305,], 2, pval_SR),
                 apply(PSS06[406:505,], 2, pval_SR))
PSS04_p <- rbind(apply(PSS04[2:101,], 2, pval_SR),
                 apply(PSS04[204:305,], 2, pval_SR),
                 apply(PSS04[406:505,], 2, pval_SR))
POW_PSS <- cbind(apply(PSS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSS04_p, 1, function(x) {mean(x<=0.05)}))


PT5S1_p <- rbind(apply(PT5S1[2:1002,], 2, pval),
                 apply(PT5S1[2006:3006,], 2, pval),
                 apply(PT5S1[4010:5010,], 2, pval))
PT5S08_p <- rbind(apply(PT5S08[2:1002,], 2, pval),
                  apply(PT5S08[2006:3006,], 2, pval),
                  apply(PT5S08[4010:5010,], 2, pval))
PT5S06_p <- rbind(apply(PT5S06[2:1002,], 2, pval),
                  apply(PT5S06[2006:3006,], 2, pval),
                  apply(PT5S06[4010:5010,], 2, pval))
PT5S04_p <- rbind(apply(PT5S04[2:1002,], 2, pval),
                  apply(PT5S04[2006:3006,], 2, pval),
                  apply(PT5S04[4010:5010,], 2, pval))
POW_PT5S <- cbind(apply(PT5S1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5S08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5S06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5S04_p, 1, function(x) {mean(x<=0.05)}))



PT4S1_p <- rbind(apply(PT4S1[2:1002,], 2, pval),
                 apply(PT4S1[2006:3006,], 2, pval),
                 apply(PT4S1[4010:5010,], 2, pval))
PT4S08_p <- rbind(apply(PT4S08[2:1002,], 2, pval),
                  apply(PT4S08[2006:3006,], 2, pval),
                  apply(PT4S08[4010:5010,], 2, pval))
PT4S06_p <- rbind(apply(PT4S06[2:1002,], 2, pval),
                  apply(PT4S06[2006:3006,], 2, pval),
                  apply(PT4S06[4010:5010,], 2, pval))
PT4S04_p <- rbind(apply(PT4S04[2:1002,], 2, pval),
                  apply(PT4S04[2006:3006,], 2, pval),
                  apply(PT4S04[4010:5010,], 2, pval))
POW_PT4S <- cbind(apply(PT4S1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4S08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4S06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4S04_p, 1, function(x) {mean(x<=0.05)}))




PT3S1_p <- rbind(apply(PT3S1[2:1002,], 2, pval),
                 apply(PT3S1[2006:3006,], 2, pval),
                 apply(PT3S1[4010:5010,], 2, pval))
PT3S08_p <- rbind(apply(PT3S08[2:1002,], 2, pval),
                  apply(PT3S08[2006:3006,], 2, pval),
                  apply(PT3S08[4010:5010,], 2, pval))
PT3S06_p <- rbind(apply(PT3S06[2:1002,], 2, pval),
                  apply(PT3S06[2006:3006,], 2, pval),
                  apply(PT3S06[4010:5010,], 2, pval))
PT3S04_p <- rbind(apply(PT3S04[2:1002,], 2, pval),
                  apply(PT3S04[2006:3006,], 2, pval),
                  apply(PT3S04[4010:5010,], 2, pval))
POW_PT3S <- cbind(apply(PT3S1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT3S08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT3S06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT3S04_p, 1, function(x) {mean(x<=0.05)}))




PT2S1_p <- rbind(apply(PT2S1[2:1002,], 2, pval),
                 apply(PT2S1[2006:3006,], 2, pval),
                 apply(PT2S1[4010:5010,], 2, pval))
PT2S08_p <- rbind(apply(PT2S08[2:1002,], 2, pval),
                  apply(PT2S08[2006:3006,], 2, pval),
                  apply(PT2S08[4010:5010,], 2, pval))
PT2S06_p <- rbind(apply(PT2S06[2:1002,], 2, pval),
                  apply(PT2S06[2006:3006,], 2, pval),
                  apply(PT2S06[4010:5010,], 2, pval))
PT2S04_p <- rbind(apply(PT2S04[2:1002,], 2, pval),
                  apply(PT2S04[2006:3006,], 2, pval),
                  apply(PT2S04[4010:5010,], 2, pval))
POW_PT2S <- cbind(apply(PT2S1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT2S08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT2S06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT2S04_p, 1, function(x) {mean(x<=0.05)}))



PPB1_p <- rbind(apply(PPB1[2:1002,], 2, pval),
                apply(PPB1[2006:3006,], 2, pval),
                apply(PPB1[4010:5010,], 2, pval))
PPB08_p <- rbind(apply(PPB08[2:1002,], 2, pval),
                 apply(PPB08[2006:3006,], 2, pval),
                 apply(PPB08[4010:5010,], 2, pval))
PPB06_p <- rbind(apply(PPB06[2:1002,], 2, pval),
                 apply(PPB06[2006:3006,], 2, pval),
                 apply(PPB06[4010:5010,], 2, pval))
PPB04_p <- rbind(apply(PPB04[2:1002,], 2, pval),
                 apply(PPB04[2006:3006,], 2, pval),
                 apply(PPB04[4010:5010,], 2, pval))
POW_PPB <- cbind(apply(PPB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PPB04_p, 1, function(x) {mean(x<=0.05)}))




PTB1_p <- rbind(apply(PTB1[2:1002,], 2, pval),
                apply(PTB1[2006:3006,], 2, pval),
                apply(PTB1[4010:5010,], 2, pval))
PTB08_p <- rbind(apply(PTB08[2:1002,], 2, pval),
                 apply(PTB08[2006:3006,], 2, pval),
                 apply(PTB08[4010:5010,], 2, pval))
PTB06_p <- rbind(apply(PTB06[2:1002,], 2, pval),
                 apply(PTB06[2006:3006,], 2, pval),
                 apply(PTB06[4010:5010,], 2, pval))
PTB04_p <- rbind(apply(PTB04[2:1002,], 2, pval),
                 apply(PTB04[2006:3006,], 2, pval),
                 apply(PTB04[4010:5010,], 2, pval))
POW_PTB <- cbind(apply(PTB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PTB04_p, 1, function(x) {mean(x<=0.05)}))




PSB1_p <- rbind(apply(PSB1[2:101,], 2, pval_SR),
                apply(PSB1[204:305,], 2, pval_SR),
                apply(PSB1[406:505,], 2, pval_SR))
PSB08_p <- rbind(apply(PSB08[2:101,], 2, pval_SR),
                 apply(PSB08[204:305,], 2, pval_SR),
                 apply(PSB08[406:505,], 2, pval_SR))
PSB06_p <- rbind(apply(PSB06[2:101,], 2, pval_SR),
                 apply(PSB06[204:305,], 2, pval_SR),
                 apply(PSB06[406:505,], 2, pval_SR))
PSB04_p <- rbind(apply(PSB04[2:101,], 2, pval_SR),
                 apply(PSB04[204:305,], 2, pval_SR),
                 apply(PSB04[406:505,], 2, pval_SR))
POW_PSB <- cbind(apply(PSB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(PSB04_p, 1, function(x) {mean(x<=0.05)}))




PT8B1_p <- rbind(apply(PT8B1[2:1002,], 2, pval),
                 apply(PT8B1[2006:3006,], 2, pval),
                 apply(PT8B1[4010:5010,], 2, pval))
PT8B08_p <- rbind(apply(PT8B08[2:1002,], 2, pval),
                  apply(PT8B08[2006:3006,], 2, pval),
                  apply(PT8B08[4010:5010,], 2, pval))
PT8B06_p <- rbind(apply(PT8B06[2:1002,], 2, pval),
                  apply(PT8B06[2006:3006,], 2, pval),
                  apply(PT8B06[4010:5010,], 2, pval))
PT8B04_p <- rbind(apply(PT8B04[2:1002,], 2, pval),
                  apply(PT8B04[2006:3006,], 2, pval),
                  apply(PT8B04[4010:5010,], 2, pval))
POW_PT8B <- cbind(apply(PT8B1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT8B08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT8B06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT8B04_p, 1, function(x) {mean(x<=0.05)}))



PT6B1_p <- rbind(apply(PT6B1[2:1002,], 2, pval),
                 apply(PT6B1[2006:3006,], 2, pval),
                 apply(PT6B1[4010:5010,], 2, pval))
PT6B08_p <- rbind(apply(PT6B08[2:1002,], 2, pval),
                  apply(PT6B08[2006:3006,], 2, pval),
                  apply(PT6B08[4010:5010,], 2, pval))
PT6B06_p <- rbind(apply(PT6B06[2:1002,], 2, pval),
                  apply(PT6B06[2006:3006,], 2, pval),
                  apply(PT6B06[4010:5010,], 2, pval))
PT6B04_p <- rbind(apply(PT6B04[2:1002,], 2, pval),
                  apply(PT6B04[2006:3006,], 2, pval),
                  apply(PT6B04[4010:5010,], 2, pval))
POW_PT6B <- cbind(apply(PT6B1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT6B08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT6B06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT6B04_p, 1, function(x) {mean(x<=0.05)}))




PT5B1_p <- rbind(apply(PT5B1[2:1002,], 2, pval),
                 apply(PT5B1[2006:3006,], 2, pval),
                 apply(PT5B1[4010:5010,], 2, pval))
PT5B08_p <- rbind(apply(PT5B08[2:1002,], 2, pval),
                  apply(PT5B08[2006:3006,], 2, pval),
                  apply(PT5B08[4010:5010,], 2, pval))
PT5B06_p <- rbind(apply(PT5B06[2:1002,], 2, pval),
                  apply(PT5B06[2006:3006,], 2, pval),
                  apply(PT5B06[4010:5010,], 2, pval))
PT5B04_p <- rbind(apply(PT5B04[2:1002,], 2, pval),
                  apply(PT5B04[2006:3006,], 2, pval),
                  apply(PT5B04[4010:5010,], 2, pval))
POW_PT5B <- cbind(apply(PT5B1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5B08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5B06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT5B04_p, 1, function(x) {mean(x<=0.05)}))




PT4B1_p <- rbind(apply(PT4B1[2:1002,], 2, pval),
                 apply(PT4B1[2006:3006,], 2, pval),
                 apply(PT4B1[4010:5010,], 2, pval))
PT4B08_p <- rbind(apply(PT4B08[2:1002,], 2, pval),
                  apply(PT4B08[2006:3006,], 2, pval),
                  apply(PT4B08[4010:5010,], 2, pval))
PT4B06_p <- rbind(apply(PT4B06[2:1002,], 2, pval),
                  apply(PT4B06[2006:3006,], 2, pval),
                  apply(PT4B06[4010:5010,], 2, pval))
PT4B04_p <- rbind(apply(PT4B04[2:1002,], 2, pval),
                  apply(PT4B04[2006:3006,], 2, pval),
                  apply(PT4B04[4010:5010,], 2, pval))
POW_PT4B <- cbind(apply(PT4B1_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4B08_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4B06_p, 1, function(x) {mean(x<=0.05)}),
                  apply(PT4B04_p, 1, function(x) {mean(x<=0.05)}))




POB1_p <- rbind(apply(POB1[2:1002,], 2, pval),
                apply(POB1[2006:3006,], 2, pval),
                apply(POB1[4010:5010,], 2, pval))
POB08_p <- rbind(apply(POB08[2:1002,], 2, pval),
                 apply(POB08[2006:3006,], 2, pval),
                 apply(POB08[4010:5010,], 2, pval))
POB06_p <- rbind(apply(POB06[2:1002,], 2, pval),
                 apply(POB06[2006:3006,], 2, pval),
                 apply(POB06[4010:5010,], 2, pval))
POB04_p <- rbind(apply(POB04[2:1002,], 2, pval),
                 apply(POB04[2006:3006,], 2, pval),
                 apply(POB04[4010:5010,], 2, pval))
POW_POB <- cbind(apply(POB1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POB08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POB06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POB04_p, 1, function(x) {mean(x<=0.05)}))


POS1_p <- rbind(apply(POS1[2:1002,], 2, pval),
                apply(POS1[2006:3006,], 2, pval),
                apply(POS1[4010:5010,], 2, pval))
POS08_p <- rbind(apply(POS08[2:1002,], 2, pval),
                 apply(POS08[2006:3006,], 2, pval),
                 apply(POS08[4010:5010,], 2, pval))
POS06_p <- rbind(apply(POS06[2:1002,], 2, pval),
                 apply(POS06[2006:3006,], 2, pval),
                 apply(POS06[4010:5010,], 2, pval))
POS04_p <- rbind(apply(POS04[2:1002,], 2, pval),
                 apply(POS04[2006:3006,], 2, pval),
                 apply(POS04[4010:5010,], 2, pval))
POW_POS <- cbind(apply(POS1_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POS08_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POS06_p, 1, function(x) {mean(x<=0.05)}),
                 apply(POS04_p, 1, function(x) {mean(x<=0.05)}))


POW_PS <- rbind(POW_POS, POW_PLS,POW_PTS, POW_PT2S, POW_PT3S, POW_PT4S, POW_PT5S, POW_PSS)
colnames(POW_PS) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_PS) <- c("POS_P", "POS_K", "POS_T",
                      "PPS_P", "PPS_K", "PPS_T",
                      "PTS_P", "PTS_K", "PTS_T",
                      "PT2S_P", "PT2S_K", "PT2S_T",
                      "PT3S_P", "PT3S_K", "PT3S_T",
                      "PT4S_P", "PT4S_K", "PT4S_T",
                      "PT5S_P", "PT5S_K", "PT5S_T",
                      "PSS_P", "PSS_K", "PSS_T")

write.csv(POW_PS, "POW_PS.csv")


POW_PB <- rbind(POW_POB, POW_PPB,POW_PTB, POW_PT4B, POW_PT5B, POW_PT6B, POW_PT8B, POW_PSB)
colnames(POW_PB) <- c("1", "0.8", "0.6", "0.4")
rownames(POW_PB) <- c("POB_P", "POB_K", "POB_T",
                      "PPB_P", "PPB_K", "PPB_T",
                      "PTB_P", "PTB_K", "PTB_T",
                      "PT4B_P", "PT4B_K", "PT4B_T",
                      "PT5B_P", "PT5B_K", "PT5B_T",
                      "PT6B_P", "PT6B_K", "PT6B_T",
                      "PT8B_P", "PT8B_K", "PT8B_T",
                      "PSB_P", "PSB_K", "PSB_T")

write.csv(POW_PB, "POW_PB.csv")
