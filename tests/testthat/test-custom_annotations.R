test_that("Custom annotation works", {

 dt_raw <- behavr::toy_ethoscope_data()
 dt <- sum_movement_detector(dt_raw, threshold=0.5)
 expect_equal(
         sum(10**((dt_raw$xy_dist_log10x1000)/1000)[1:20]),
         dt$sum_movement[1]
 )

})
