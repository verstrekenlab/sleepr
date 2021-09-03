test_that("Custom annotation works", {

 dt_raw <- behavr::toy_ethoscope_data()
 dt <- sum_movement_detector(dt_raw, threshold=0.5)
 expect_equal(
         dt_raw[1:20, sum(10**(xy_dist_log10x1000/1000))],
         dt[1, sum_movement]
 )

})
