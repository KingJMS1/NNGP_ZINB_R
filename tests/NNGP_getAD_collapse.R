source("../R/NNGP_getAD_collapse.R")
neardist = read.table("data/neardist")
neardistM = read.table("data/neardistM")
N = 30
M = 10
phi = 8.336574

output = getAD(neardist, neardistM, N, M, phi)
# print(output)
correct_output = read.table("data/getAD")
output[is.na(output)] = -999
correct_output[is.na(correct_output)] = -999
print(sum(abs(output - correct_output)))