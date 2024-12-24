# -------------------------------------------------------------------------------------- #
#                                    data preparation                                    #
# -------------------------------------------------------------------------------------- #

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
catalogue <- read.csv("16runs_catalogue.csv")

suppressPackageStartupMessages(library(dplyr))

# calculate the number of all regular designs
orig <- catalogue %>% 
  filter(num_of_columns >= 6, num_of_columns <= 7) %>%
  nrow()
comp <- catalogue %>% 
  filter(num_of_columns >= 5, num_of_columns <= 7) %>%
  nrow()
cat(paste0("Num of regular designs: ", orig + comp, "\n"))

indep_cols <- read.csv("16runs_indep.csv") |>
  as.matrix()
generator <- read.csv("16runs_generator.csv") |>
  as.matrix()
lines <- read.csv("16runs_lines.csv") |>
  as.matrix()

# load in design A of SOA 2+
# wlp in good_A starts at length = 3
good_A <- read.csv("16runs_good_A.csv")
good_A <- good_A %>%
  filter(num_of_columns >= 6)
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

is_good_A <- function(a_columns) {
  idx <- which((matrix(lines %in% a_columns, nrow(lines)) %*% rep(1,3)) == 1)
  return(all(a_columns %in% lines[idx, ]))
}

get_b_set <- function(a_columns, s) {
  possible_b <- list()
  for (i in 1:length(a_columns)) {
    idx <- which((lines == a_columns[i]) %*% rep(1, (s + 1)) == 1)
    nums_of_col_in_a <- matrix(lines[idx, ] %in% a_columns, length(idx)) %*% rep(1, (s + 1))
    exactly_one <- (lines[idx, ])[which(nums_of_col_in_a == 1), ]
    possible_b[[i]] <- setdiff(unique(as.vector(exactly_one)), a_columns[i])
  }
  return(possible_b)
}

s22 <- function(d, s){
  # count: if stratififed, c(0,0) should appear $count times
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

# -------------------------------------------------------------------------------------- #
#                                     2x2x2 then 4x4                                     #
# -------------------------------------------------------------------------------------- #

cat("----- Optimal alpha-beta designs -----\n")
# first maximize 2x2x2
filtered <- data.frame(matrix(nrow = 0, ncol = 5))
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

for (k in 1:nrow(filtered)) {

  cat(paste0("Processing s22: ", k, "/", nrow(filtered), "\r"))
  a_columns <- str_to_vec(filtered$columns[k])
  A <- (indep_cols %*% generator[, a_columns]) %% 2
  b_set <- get_b_set(a_columns, 2)
  
  m <- length(b_set)

  max_positions <- lapply(b_set, length) |> unlist()
  curr_position <- rep(1, m)
  
  b_columns <- rep(1, m)
  best_b_columns <- c()
  max_count <- 0
  
  while (TRUE) {
    # do things 
    for (j in 1:m) {
      b_columns[j] <- b_set[[j]][curr_position[j]]
    }
    B <- (indep_cols %*% generator[, b_columns]) %% 2
    D <- 2 * A + B
    count <- count_s22(D, 2)
    if (count > max_count) {
      best_b_columns <- b_columns
      max_count <- count
    }
    if (all(max_positions == curr_position)) break
    # move to the next design
    cursor <- 1
    while (curr_position[cursor] + 1 > max_positions[cursor]) {
      curr_position[cursor] <- 1
      cursor <- cursor + 1
    }
    curr_position[cursor] <- curr_position[cursor] + 1
  }
  
  good_B[k] <- vec_to_str(best_b_columns)
  s22_max[k] <- max_count
}
cat("\n")

result <- data.frame(b_columns=good_B, s22_max=s22_max)
result <- cbind(filtered, result)
cat(paste0("After s22: ", nrow(result), "\n"))

# final processing
result$s111 <- rep(0, nrow(result))
for (i in 1:nrow(result)) {
  # B can be ignored when counting s111
  A <- (indep_cols %*% generator[, str_to_vec(result$columns[i])]) %% 2
  D <- 2 * A
  result$s111[i] <- count_s111(D, 2)
}
result <- result %>%
  mutate(s111_full = choose(num_of_columns, 3), s22_full = choose(num_of_columns, 2)) %>%
  mutate(s111_prop = (s111 / s111_full), s22_prop = (s22_max / s22_full)) %>%
  select(-c(s111_full, s22_full))
result <- result[, c(1, 2, 3, 6, 8, 9, 7, 10)]
colnames(result) <- c("idx", "m", "a_columns", "b_columns", "s111", "s111_prop", "s22", "s22_prop")

write.csv(result, "16runs_alpha_beta.csv", row.names = F)

# -------------------------------------------------------------------------------------- #
#                                     4x4 then 2x2x2                                     #
# -------------------------------------------------------------------------------------- #

cat("----- Optimal beta-alpha designs -----\n")
# first maximize 4x4
good_B <- c()
s22_max <- c()

for (i in 1:nrow(good_A)) {

  cat(paste0("Processing s22: ", i, "/", nrow(good_A), "\r"))
  a_columns <- str_to_vec(good_A$columns[i])
  A <- (indep_cols %*% generator[, a_columns]) %% 2
  b_set <- get_b_set(a_columns, 2)
  
  m <- length(b_set)
  
  max_positions <- lapply(b_set, length) |> unlist()
  curr_position <- rep(1, m)
  
  b_columns <- rep(1, m)
  best_b_columns <- c()
  max_count <- 0
  
  while (TRUE) {
    # do things 
    for (j in 1:m) {
      b_columns[j] <- b_set[[j]][curr_position[j]]
    }
    B <- (indep_cols %*% generator[, b_columns]) %% 2
    D <- 2 * A + B
    count <- count_s22(D, 2)
    if (count > max_count) {
      best_b_columns <- b_columns
      max_count <- count
    }
    if (all(max_positions == curr_position)) break
    # move to the next design
    cursor <- 1
    while (curr_position[cursor] + 1 > max_positions[cursor]) {
      curr_position[cursor] <- 1
      cursor <- cursor + 1
    }
    curr_position[cursor] <- curr_position[cursor] + 1
  }
  good_B[i] <- vec_to_str(best_b_columns)
  s22_max[i] <- max_count
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
  A <- (indep_cols %*% generator[, str_to_vec(filtered$a_columns[i])]) %% 2
  D <- 2 * A
  filtered$s111[i] <- count_s111(D, 2)
}
filtered <- filtered %>%
  mutate(s111_full = choose(num_of_columns, 3), s22_full = choose(num_of_columns, 2)) %>%
  mutate(s111_prop = (s111 / s111_full), s22_prop = (s22_max / s22_full)) %>%
  select(-c(s111_full, s22_full))
filtered <- filtered[, c(1, 2, 3, 6, 8, 9, 7, 10)]
colnames(filtered) <- c("idx", "m", "a_columns", "b_columns", "s111", "s111_prop", "s22", "s22_prop")

write.csv(filtered, "16runs_beta_alpha.csv", row.names = F)
