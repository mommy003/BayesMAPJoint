
###test on a simple chromosome example
m <- nrow(LD_exact)

chr <- c(
  rep("1", 100),
  rep("2", 100)
)

position <- c(
  seq(1, 1000000, length.out = 100),
  seq(1, 1000000, length.out = 100)
)

block_table <- data.frame(
  chr = c("1", "2"),
  start = c(1, 1),
  end = c(1000000, 1000000)
)


block_info <- make_ld_blocks(
  chr = chr,
  position = position,
  block_table = block_table
)

block_indices_chr <- block_info$block_indices

LD_blocks_chr <- extract_ld_blocks(
  LD = LD_exact,
  block_indices = block_indices_chr
)

check_ld_block_quality(
  LD = LD_exact,
  block_indices = block_indices_chr
)
