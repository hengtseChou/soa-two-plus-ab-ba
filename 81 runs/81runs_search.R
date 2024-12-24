# -------------------------------------------------------------------------------------- #
#                                    data preparation                                    #
# -------------------------------------------------------------------------------------- #

set.seed(8832)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
catalogue <- read.csv("81runs_catalogue.csv")

suppressPackageStartupMessages(library(dplyr))

# calculate the number of all regular designs
orig <- catalogue %>% 
  filter(num_of_columns >= 11, num_of_columns <= 20) %>%
  nrow()
# 重複的只要做一邊就可以
comp <- catalogue %>% 
  filter(num_of_columns >= 11, num_of_columns <= 19) %>%
  nrow()
cat(paste0("Num of regular designs: ", orig + comp, "\n"))

indep_cols <- read.csv("81runs_indep.csv") |>
  as.matrix()
generator <- read.csv("81runs_generator.csv") |>
  as.matrix()
lines <- read.csv("81runs_lines.csv") |>
  as.matrix()

# load in design A of SOA 2+
# wlp in good_A starts at length = 3
good_A <- read.csv("81runs_good_A.csv") %>%
  filter(num_of_columns >= 11)
cat(paste0("Num of good A: ", nrow(good_A), "\n"))

# -------------------------------------------------------------------------------------- #
#                                 functions and algorithm                                #
# -------------------------------------------------------------------------------------- #

str_to_vec <- function(s) {
  return(as.numeric(unlist(strsplit(s, " "))))
}

vec_to_str <- function(v) {
  if (is.null(v)) return(NA)
  return(paste(v, collapse = " "))
}

get_b_set <- function(a_columns, s) {
  possible_b <- list()
  for (i in 1:length(a_columns)) {
    idx <- which((lines == a_columns[i]) %*% rep(1, (s+1)) == 1)
    nums_of_col_in_a <- matrix(lines[idx, ] %in% a_columns, length(idx)) %*% rep(1, (s+1))
    exactly_one <- (lines[idx, ])[which(nums_of_col_in_a == 1), ]
    possible_b[[i]] <- setdiff(unique(as.vector(exactly_one)), a_columns[i])
  }
  return(possible_b)
}

s22 <- function(d, s){
  # count: if stratified, c(0,0) should appear $count times
  count <- nrow(d) / s ^ 4
  pairs <- combn(ncol(d), 2)
  has_s22 <- rep(1, ncol(pairs))
  for (i in 1:ncol(pairs)) {
    pair <- d[ ,pairs[ ,i]]
    if (sum(pair %*% c(1,1) == 0) != count) has_s22[i] <- 0
  }
  return(has_s22)
}

count_s22 <- function(d, s) {
  return(sum(s22(d, s)))
}

n_dis <- function(x) { #number of distinct elements
  length(unique(x))
}

s111 <- function(d, s) {
  A <- combn(ncol(d), 3)
  check <- rep(0, ncol(A))
  for(i in 1:ncol(A)) {
    dp <- d[,A[,i]]
    if (n_dis(floor(dp / s) %*% s^(0:2)) == s^3) check[i] <- 1
  }
  return(check)
}

count_s111 <- function(d, s) {
  return(sum(s111(d, s)))
}

generate_random_b_cols <- function(b_set) {
  b_columns <- c()
  for (i in 1:length(b_set)) b_columns[i] <- sample(b_set[[i]], 1)
  return(b_columns)
}

greedy <- function(a_columns, iter, s) {
  A <- (indep_cols %*% generator[, a_columns]) %% s
  b_set <- get_b_set(a_columns, s)
  b_columns <- generate_random_b_cols(b_set)
  for (t in 1:iter) {
    for (i in 1:length(b_set)) {
      counts <- c()
      for (j in 1:length(b_set[[i]])) {
        tmp_b_cols <- b_columns
        tmp_b_cols[i] <- (b_set[[i]])[j]
        tmp_B <- (indep_cols %*% generator[, tmp_b_cols]) %% s
        tmp_D <- s * A + tmp_B
        counts[j] <- count_s22(tmp_D, s)
      }
      b_columns[i] <- (b_set[[i]])[which.max(counts)]
    }
  }
  return(b_columns)
}

# -------------------------------------------------------------------------------------- #
#                                     3x3x3 then 9x9                                     #
# -------------------------------------------------------------------------------------- #

cat("----- Optimal alpha-beta designs -----\n")
# first maximize 2x2x2
filtered <- data.frame(matrix(nrow = 0, ncol = 5))
library(dplyr)
# original designs: filter by the least words with length 3
range1 <- good_A %>% 
  filter(is_comp == F) %>% 
  select(num_of_columns) %>% 
  range
for (m in seq(range1[1], range1[2])) {
  min_wlp <- 10000
  entries <- good_A %>% filter(is_comp == F, num_of_columns == m)
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] > min_wlp) next
    if (wlp[1] < min_wlp) min_wlp <- wlp[1]
  }
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] == min_wlp) filtered <- rbind(filtered, entries[i, ])
  }
}
# comp. design: filter by the most words with length 3
range2 <- good_A %>% 
  filter(is_comp == T) %>% 
  select(num_of_columns) %>% 
  range
for (m in seq(range2[1], range2[2])) {
  max_wlp <- 0
  entries <- good_A %>% 
    filter(is_comp == T, num_of_columns == m) %>% 
    arrange(desc(row_number()))
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] < max_wlp) next
    if (wlp[1] > max_wlp) max_wlp <- wlp[1]
  }
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] == max_wlp) filtered <- rbind(filtered, entries[i, ])
  }
}
cat(paste0("After s111: ", nrow(filtered), "\n"))

# then maximize 4x4
good_B <- c()
s22_max <- c()
iter <- 5

for (k in 1:nrow(filtered)) {

  cat(paste0("Processing s22: ", k, "/", nrow(filtered), "\r"))
  a_columns <- str_to_vec(filtered$columns[k])
  A <- (indep_cols %*% generator[, a_columns]) %% 3
  b_columns <- greedy(a_columns, iter, 3)
  B <- (indep_cols %*% generator[, b_columns]) %% 3
  D <- 3 * A + B
  
  good_B[k] <- vec_to_str(b_columns)
  s22_max[k] <- count_s22(D, 3)
}
cat("\n")

result <- data.frame(b_columns=good_B, s22_max=s22_max)
result <- cbind(filtered, result)
result <- result %>%
  group_by(num_of_columns) %>%
  filter(s22_max == max(s22_max, na.rm=TRUE)) %>%
  mutate(b_columns = ifelse(s22_max == 0, NA, b_columns))
cat(paste0("After s22: ", nrow(result), "\n"))

# final processing
result$s111 <- rep(0, nrow(result))
for (i in 1:nrow(result)) {
  # B can be ignored when counting s111
  A <- (indep_cols %*% generator[, str_to_vec(result$columns[i])]) %% 3
  D <- 3 * A
  result$s111[i] <- count_s111(D, 3)
}
result <- result %>%
  mutate(s111_full = choose(num_of_columns, 3), s22_full = choose(num_of_columns, 2)) %>%
  mutate(s111_prop = (s111 / s111_full), s22_prop = (s22_max / s22_full)) %>%
  select(-c(s111_full, s22_full))
result <- result[, c(1, 2, 3, 6, 8, 9, 7, 10)]
colnames(result) <- c("idx", "m", "a_columns", "b_columns", "s111", "s111_prop", "s22", "s22_prop")

write.csv(result, "81runs_alpha_beta.csv", row.names = F)

# -------------------------------------------------------------------------------------- #
#                                     9x9 then 3x3x3                                     #
# -------------------------------------------------------------------------------------- #

cat("----- Optimal beta-alpha designs -----\n")
# first maximize 4x4
good_B <- c()
s22_max <- c()
iter <- 5

for (k in 1:nrow(good_A)) {

  cat(paste0("Processing s22: ", k, "/", nrow(good_A), "\r"))
  a_columns <- str_to_vec(good_A$columns[k])
  A <- (indep_cols %*% generator[, a_columns]) %% 3
  b_columns <- greedy(a_columns, iter, 3)
  B <- (indep_cols %*% generator[, b_columns]) %% 3
  D <- 3 * A + B
  
  good_B[k] <- vec_to_str(b_columns)
  s22_max[k] <- count_s22(D, 3)
}
cat("\n")

result <- data.frame(b_columns=good_B, s22_max=s22_max)
result <- cbind(good_A, result)
colnames(result)[3] <- "a_columns"

result_s22_max <- result %>%
  group_by(num_of_columns) %>%
  filter(s22_max == max(s22_max, na.rm=TRUE))
cat(paste0("After s22: ", nrow(result_s22_max), "\n"))

# then 2x2x2
filtered <- data.frame(matrix(nrow = 0, ncol = 7))
# original designs: filter by the least words with length 3
range1 <- result_s22_max %>% 
  filter(is_comp == F) %>% 
  select(num_of_columns) %>% 
  range
for (m in seq(range1[1], range1[2])) {
  min_wlp <- 10000
  entries <- result_s22_max %>% filter(is_comp == F, num_of_columns == m)
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] > min_wlp) next
    if (wlp[1] < min_wlp) min_wlp <- wlp[1]
  }
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] == min_wlp) filtered <- rbind(filtered, entries[i, ])
  }
}
# comp. design: filter by the most words with length 3
range2 <- result_s22_max %>% 
  filter(is_comp == T) %>% 
  select(num_of_columns) %>% 
  range
for (m in seq(range2[1], range2[2])) {
  max_wlp <- 0
  entries <- result_s22_max %>% 
    filter(is_comp == T, num_of_columns == m) %>% 
    arrange(desc(row_number()))
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] < max_wlp) next
    if (wlp[1] > max_wlp) max_wlp <- wlp[1]
  }
  for (i in 1:nrow(entries)) {
    wlp <- str_to_vec(entries$wlp[i])
    if (wlp[1] == max_wlp) filtered <- rbind(filtered, entries[i, ])
  }
}
cat(paste0("After s111: ", nrow(filtered), "\n"))

# final processing
filtered$s111 <- rep(0, nrow(filtered))
for (i in 1:nrow(filtered)) {
  # B can be ignored when counting s111
  A <- (indep_cols %*% generator[, str_to_vec(filtered$a_columns[i])]) %% 3
  D <- 3 * A
  filtered$s111[i] <- count_s111(D, 3)
}
filtered <- filtered %>%
  mutate(s111_full = choose(num_of_columns, 3), s22_full = choose(num_of_columns, 2)) %>%
  mutate(s111_prop = (s111 / s111_full), s22_prop = (s22_max / s22_full)) %>%
  select(-c(s111_full, s22_full))
filtered <- filtered[, c(1, 2, 3, 6, 8, 9, 7, 10)]
colnames(filtered) <- c("idx", "m", "a_columns", "b_columns", "s111", "s111_prop", "s22", "s22_prop")

write.csv(filtered, "81runs_beta_alpha.csv", row.names = F)
