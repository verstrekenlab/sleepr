test_that("Custom annotation works", {

 dt <- behavr::toy_ethoscope_data()
 distance_sum_enclosed(dt)

 # dt1 <- sleepr::custom_annotation_wrapper(sleepr::distance_sum_enclosed)(dt)
 dt2 <- sum_movement_detector(dt, threshold=0.5)

 sleep_annotation(dt)
 args <- list(data=dt, FUN=list(sum_movement_detector), time_window_length=10)
 dt <- do.call(scopr::annotate_all, args)
 expect_true("micromovement" %in% colnames(dt))

})
