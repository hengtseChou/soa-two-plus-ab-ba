# -------------------------------------------------------------------------------------- #
#                                       indep cols                                       #
# -------------------------------------------------------------------------------------- #

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
e1 <- rep(0:1, each=16)
e2 <- rep(0:1, each=8, 2)
e3 <- rep(0:1, each=4, 4)
e4 <- rep(0:1, each=2, 8)
e5 <- rep(0:1, each=1, 16)
indep_cols <- cbind(e1, e2, e3, e4, e5)
colnames(indep_cols) <- 1:5
write.csv(indep_cols, "32runs_indep.csv", row.names=F)

# -------------------------------------------------------------------------------------- #
#                                        generator                                       #
# -------------------------------------------------------------------------------------- #

full <- function(l, n) { # return an l^n full factorial design
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
generator <- t(full(2, 5)[-1, 5:1])
write.csv(generator, "32runs_generator.csv", row.names=F)

# -------------------------------------------------------------------------------------- #
#                                          lines                                         #
# -------------------------------------------------------------------------------------- #

lines <- c()
C2 <- combn(31, 2)
for (i in 1:ncol(C2)) {
  s2 <- (generator[, C2[1, i]] + generator[, C2[2, i]]) %% 2
  tmp <- which(2 ^ (0:4) %*% generator == sum(s2 * 2 ^ (0:4)))
  lines <- rbind(lines, c(C2[, i], tmp))
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
write.csv(unique_lines, "32runs_lines.csv", row.names=F)

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
  idx <- which((matrix(lines %in% a_columns, nrow(lines)) %*% rep(1,3)) == 1)
  return(all(a_columns %in% lines[idx, ]))
}

saturated <- 1:31
catalogue <- read.csv("32runs_catalogue.csv")
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
  if (is_good_A(complementary)) {
    comp_good_A_idx <- c(comp_good_A_idx, i)
  }
}

library(stringr)
orig_good_A <- c()
orig_good_A_num_of_columns <- c()
orig_good_A_wlp <- c()
for (idx in orig_good_A_idx) {
  orig_good_A <- c(orig_good_A, catalogue$columns[idx])
  orig_good_A_num_of_columns <- c(orig_good_A_num_of_columns, catalogue$num_of_columns[idx])
  orig_good_A_wlp <- c(orig_good_A_wlp, str_sub(catalogue$wlp[idx], 5))
}

comp_good_A <- c()
comp_good_A_num_of_columns <- c()
comp_good_A_wlp <- c()
for (idx in comp_good_A_idx) {
  original <- str_to_vec(catalogue$columns[idx])
  comp <- saturated[-original]
  comp_good_A_num_of_columns <- c(comp_good_A_num_of_columns, length(comp))
  comp_good_A <- c(comp_good_A, vec_to_str(comp))
  comp_good_A_wlp <- c(comp_good_A_wlp, str_sub(catalogue$wlp[idx], 5))
}

good_A <- data.frame(idx=c(orig_good_A_idx, comp_good_A_idx), 
                     num_of_columns=c(orig_good_A_num_of_columns, comp_good_A_num_of_columns), 
                     columns=c(orig_good_A, comp_good_A),
                     wlp=c(orig_good_A_wlp, comp_good_A_wlp),
                     is_comp=c(rep(F, length(orig_good_A)), rep(T, length(comp_good_A))))

library(dplyr)
good_A <- good_A %>% arrange(num_of_columns, idx)
write.csv(good_A, "32runs_good_A.csv", row.names = FALSE)
