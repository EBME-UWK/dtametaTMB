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

test_that("RevMan sample-size transformation is reproduced", {
  dat <- data.frame(
    n1 = c(25, 100, 400),
    n0 = c(25, 100, 400)
  )
  
  out <- getWEIGHTS(dat, "sampsize_revman")
  
  expect_equal(
    out$se * 5,
    c(7.2, 10.4, 16.8)
  )
  
  expect_equal(
    out$sp * 5,
    c(7.2, 10.4, 16.8)
  )
})

test_that("RevMan inverse-SE scaling retains baseline size", {
  dat <- data.frame(
    n1 = 100,
    n0 = 100,
    sens = 1,
    spec = 1
  )
  
  out <- getWEIGHTS(dat, "se_revman")
  
  expect_equal(out$se * 5, 4)
  expect_equal(out$sp * 5, 4)
})