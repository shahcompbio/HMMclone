state <- c(replicate(10, 2), replicate(10, 4), replicate(10, 2), replicate(10, 12))
copy <- lapply(state, function(x) rnorm(1, mean = x, sd = 0.2)) %>% unlist()
df <- data.frame(clone_id = "A",
                 chr = "1",
                 start = 1:length(state),
                 end = 1:length(state) + 1,
                 copy = copy)
df_fit <- HMMclone(df, sd_value = 0.2)


test_that("Inference is correct", {
  expect_equal(df_fit$state, state)
})
