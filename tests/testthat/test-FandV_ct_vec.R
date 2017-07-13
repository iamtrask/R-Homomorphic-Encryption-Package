context("FandV scheme cipher text vectors")

test_that("Vectors", {
  p <- pars("FandV")
  keys <- keygen(p)
  ct1 <- enc(keys$pk, 2)
  ct2 <- enc(keys$pk, 3)
  ct3 <- enc(keys$pk, -4)
  
  a <- c(ct1, ct2)
  a <- c(a, ct3)
  
  b <- enc(keys$pk, c(2,3,-4))
  
  expect_that(dec(keys$sk, a[2]), equals(3))
  expect_that(dec(keys$sk, a[3:1][3]), equals(2))
  expect_that(dec(keys$sk, a[-c(2,1)]), equals(-4))
  expect_that(dec(keys$sk, b), equals(c(2,3,-4)))
  b[1] <- a[3]
  expect_that(dec(keys$sk, b), equals(c(-4,3,-4)))
  expect_that(length(a), equals(3))
})

test_that("Vector operations", {
  p <- pars("FandV")
  keys <- keygen(p)
  ct1 <- enc(keys$pk, 2)
  ct2 <- enc(keys$pk, 3)
  ct3 <- enc(keys$pk, -4)
  
  a <- c(ct1, ct2)
  a <- c(a, ct3)
  
  ct1 <- enc(keys$pk, 5)
  ct2 <- enc(keys$pk, -2)
  ct3 <- enc(keys$pk, 6)
  
  b <- c(ct2, ct3)
  b <- c(ct1, b)
  
  ct <- enc(keys$pk, 2:4)
  
  expect_that(dec(keys$sk, (a+b)[1]), equals(7))
  expect_that(dec(keys$sk, (a+b)[2]), equals(1))
  expect_that(dec(keys$sk, (a+b)[3]), equals(2))
  expect_that(dec(keys$sk, a-b), equals(c(-3, 5, -10)))
  expect_that(dec(keys$sk, c(a, ct1)-c(ct2, ct3)), equals(c(4, -3, -2, -1)))
  ### Why does including this test (which passes) result in an error if the inner product test is done after, not before?  Error in signalCondition(e): no function to return from, jumping to top level
  ### related? http://tolstoy.newcastle.edu.au/R/e16/devel/11/12/0384.html
  expect_that(dec(keys$sk, c(ct2, ct3)-c(a, ct1)), equals(c(-4, 3, 2, 1)))
  expect_that(dec(keys$sk, (a*b)[1]), equals(10))
  expect_that(dec(keys$sk, (a*b)[2]), equals(-6))
  expect_that(dec(keys$sk, (a*b)[3]), equals(-24))
  expect_that(dec(keys$sk, (a*ct1)[1]), equals(10))
  expect_that(dec(keys$sk, (a*ct1)[2]), equals(15))
  expect_that(dec(keys$sk, (a*ct1)[3]), equals(-20))
  expect_that(dec(keys$sk, sum(ct)), equals(9))
  expect_that(dec(keys$sk, prod(ct)), equals(24))
  expect_that(dec(keys$sk, a%*%b), equals(-20))
})

test_that("Vector ops (imbalanced)", {
  p <- pars("FandV")
  k <- keygen(p)
  x <- 1:10
  
  y <- x
  expect_warning(y <- y * y[5:10])
  xct <- enc(k$pk, x)
  expect_warning(xct <- xct * xct[5:10])
  expect_that(dec(k$sk, xct), equals(y))
  
  y <- x
  expect_warning(y <- y[5:10] * y)
  xct <- enc(k$pk, x)
  expect_warning(xct <- xct[5:10] * xct)
  expect_that(dec(k$sk, xct), equals(y))
})

test_that("Vector rep", {
  p <- pars("FandV")
  keys <- keygen(p)
  
  x <- 1:4
  ctx <- enc(keys$pk, x)
  
  expect_that(dec(keys$sk, rep(ctx, 2)), equals(rep(x, 2)))
  expect_that(dec(keys$sk, rep(ctx, each=2)), equals(rep(x, each=2)))
  expect_that(dec(keys$sk, rep(ctx, c(2,2,2,2))), equals(rep(x, c(2,2,2,2))))
  expect_that(dec(keys$sk, rep(ctx, c(2,1,2,1))), equals(rep(x, c(2,1,2,1))))
  expect_that(dec(keys$sk, rep(ctx, each=2, len=4)), equals(rep(x, each=2, len=4)))
  expect_that(dec(keys$sk, rep(ctx, each=2, len=10)), equals(rep(x, each=2, len=10)))
  expect_that(dec(keys$sk, rep(ctx, each=2, times=3)), equals(rep(x, each=2, times=3)))
})

test_that("Vector assignment", {
  p <- pars("FandV")
  k <- keygen(p)
  x <- 1:10
  
  y <- x
  expect_warning(y[1:3] <- y[5:10])
  xct <- enc(k$pk, x)
  expect_warning(xct[1:3] <- xct[5:10])
  expect_that(dec(k$sk, xct), equals(y))
  
  y <- x
  expect_warning(y[1:3] <- y[6:7])
  xct <- enc(k$pk, x)
  expect_warning(xct[1:3] <- xct[6:7])
  expect_that(dec(k$sk, xct), equals(y))
})
