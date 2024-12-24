# -------------------------------------------------------------------------------------- #
#                                       indep cols                                       #
# -------------------------------------------------------------------------------------- #

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
e1 <- rep(0:2, each=27)
e2 <- rep(0:2, each=9, 3)
e3 <- rep(0:2, each=3, 9)
e4 <- rep(0:2, each=1)
indep_cols <- cbind(e1, e2, e3, e4)
colnames(indep_cols) <- 1:4
write.csv(indep_cols, "81runs_indep.csv", row.names=F)

# -------------------------------------------------------------------------------------- #
#                                        generator                                       #
# -------------------------------------------------------------------------------------- #

full <- function(l, n){ # return an l^n full factorial design
  if (n == 1) return(cbind(1:l-1))
  a <- as.integer(gl(l, l^n/l, l^n))
  if (n > 1) {
    for (i in 2:n) {
      a <- cbind(a, as.integer(gl(l, l^(n-i), l^n)))
    }
  }
  d <- a - 1
  colnames(d) <- 1:ncol(d)
  return(d)
}
generator <- t(full(3, 4)[c(2, 4, 5, 8, 10, 11, 13, 14, 17, 20, 
                      22, 23, 26, 28, 29, 31, 32, 35, 37, 38,
                      40, 41, 44, 47, 49, 50, 53, 56, 58, 59, 
                      62, 64, 65, 67, 68, 71, 74, 76, 77, 80), 4:1])
write.csv(generator, "81runs_generator.csv", row.names=F)

# -------------------------------------------------------------------------------------- #
#                                          lines                                         #
# -------------------------------------------------------------------------------------- #

first_is_one <- function(x) { # check whether the first non-zero element of x is one
  for (i in 1:length(x)) {
    if (x[i] == 0) next
    if (x[i] == 1) return(T)
    if (x[i] != 1) return(F)
  }
}
mapping <- 1:40
names(mapping) <- t(generator) %*% 3^(0:3)

lines <- matrix(0, choose(40, 2), 4)
tmp <- matrix(0, 4 ,4)
k <- 1
for (i in 1:39) {
  for (j in (i + 1):40){
    tmp[1, ] <- generator[, i]
    tmp[2, ] <- generator[, j]
    tmp[3, ] <- (generator[, i] + generator[, j]) %% 3
    tmp[4, ] <- (generator[, i] + 2 * generator[, j]) %% 3
    if (!first_is_one(tmp[3, ])) tmp[3, ] <- (tmp[3, ] * 2) %% 3
    if (!first_is_one(tmp[4, ])) tmp[4, ] <- (tmp[4, ] * 2) %% 3
    lines[k, ] <- sapply(tmp %*% 3 ^ (0:3), function(x) mapping[[as.character(x)]])
    k <- k + 1
  }
}

Unique <- function(X){
  X.u <- array(0, dim=dim(X))
  X.u[1, ] <- X[1, ]
  k <- 1
  for (i in 2:nrow(X)) {
    IN <- F
    for(j in 1:k){
      if (setequal(X[i, ], X.u[j, ])) {IN <- T; break}
    }
    if (!IN) {k <- k + 1; X.u[k, ] <- X[i, ]}
  }
  return(X.u[1:k, ])
}
unique_lines <- Unique(lines)
write.csv(unique_lines, "81runs_lines.csv", row.names=F)

# -------------------------------------------------------------------------------------- #
#                                     generate good A                                    #
# -------------------------------------------------------------------------------------- #

str_to_vec <- function(s) {
  return(as.numeric(unlist(strsplit(s, " "))))
}

vec_to_str <- function(v) {
  return(paste(v, collapse = " "))
}

is_good_A <- function(a_columns) {
  idx <- which((matrix(lines %in% a_columns, nrow(lines)) %*% rep(1,4)) == 1)
  return(all(a_columns %in% lines[idx, ]))
}

saturated <- 1:40
catalogue <- read.csv("81runs_catalogue.csv")
orig_good_A_idx <- c()
comp_good_A_idx <- c()
for (i in 1:length(catalogue$columns)) {
  a_columns <- str_to_vec(catalogue$columns[i])
  complementary <- saturated[-a_columns]
  # original design
  if (is_good_A(a_columns)) {
    orig_good_A_idx <- c(orig_good_A_idx, i)
  }
  # complementary design
  if (catalogue$num_of_columns[i] == 20) {
    next
  }
  else {
    if (is_good_A(complementary)) {
      comp_good_A_idx <- c(comp_good_A_idx, i)
    }
  }
}

orig_good_A <- c()
orig_good_A_num_of_columns <- c()
orig_good_A_wlp <- c()
for (idx in orig_good_A_idx) {
  orig_good_A <- c(orig_good_A, catalogue$columns[idx])
  orig_good_A_num_of_columns <- c(orig_good_A_num_of_columns, catalogue$num_of_columns[idx])
  orig_good_A_wlp <- c(orig_good_A_wlp, catalogue$wlp[idx])
}

comp_good_A <- c()
comp_good_A_num_of_columns <- c()
comp_good_A_wlp <- c()
for (idx in comp_good_A_idx) {
  original <- str_to_vec(catalogue$columns[idx])
  comp <- saturated[-original]
  comp_good_A_num_of_columns <- c(comp_good_A_num_of_columns, length(comp))
  comp_good_A <- c(comp_good_A, vec_to_str(comp))
  comp_good_A_wlp <- c(comp_good_A_wlp, catalogue$wlp[idx])
}

good_A <- data.frame(idx=c(orig_good_A_idx, comp_good_A_idx), 
                     num_of_columns=c(orig_good_A_num_of_columns, comp_good_A_num_of_columns), 
                     columns=c(orig_good_A, comp_good_A),
                     wlp=c(orig_good_A_wlp, comp_good_A_wlp),
                     is_comp=c(rep(F, length(orig_good_A)), rep(T, length(comp_good_A))))

library(dplyr)
good_A <- good_A %>% arrange(num_of_columns, idx)
write.csv(good_A, "81runs_good_A.csv", row.names = FALSE)

