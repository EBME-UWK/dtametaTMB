test_that("getWEIGHTS supports all declared methods", {
  dat <- getXP(anticcp)
  
  methods <- c(
    "equal",
    "sampsize",
    "se",
    "sampsize_revman",
    "se_revman"
  )
  
  for (method in methods) {
    out <- getWEIGHTS(dat, method)
    
    expect_s3_class(out, "data.frame")
    expect_named(out, c("sp", "se"))
    expect_length(out$sp, nrow(dat))
    expect_length(out$se, nrow(dat))
    expect_true(all(is.finite(out$sp)))
    expect_true(all(is.finite(out$se)))
  }
})
