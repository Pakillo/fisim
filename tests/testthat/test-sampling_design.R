test_that("xy_sample", {
  #--  xy_sample_random  ------------------------------------------------------
  set.seed(2020)

  # hberg_beech
  expect_equal(
    fisim:::xy_sample_random(hberg_beech$boundary, n = 10, M = 1),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point = 1:10,
      x_s = c(79.31, 8.58, 30.31, 67.25, 66.66, 22.5, 2.81, 35.96, 36.53, 96.2),
      y_s = c(12.54, 33.49, 18.87, 116.49, 112.84, 48.87, 71.51, 103.79, 80.58,
              65.66)),
    tolerance = 5e-2)

  # kalimanton_peat
  expect_equal(
    fisim:::xy_sample_random(kalimantan_peat$boundary, n = 10, M = 1),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point = 1:10,
      x_s = c(30.31, 50.33, 21.67, 43.7, 41.25, 2.4, 42.54, 80.85, 49.79, 70.36),
      y_s = c(1.41, 67.11, 16.16, 10.41, 113.57, 5.79, 83.05, 17.72, 18.63,
              87.13)),
    tolerance = 5e-2)

  # random_tree
  expect_equal(
    fisim:::xy_sample_random(random_tree$boundary, n = 10, M = 1),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point = 1:10,
      x_s = c(348.85, 723.47, 339.19, 332.43, 672.37, 877.51, 860.89, 965.69,
              767.98, 133.71),
      y_s = c(452.27, 18.07, 276.93, 313.79, 145.1, 146.16, 172.26, 283.83,
              56.46, 236.48)),
    tolerance = 5e-2)


  #--  xy_sample_regular  -----------------------------------------------------
  set.seed(2020)

  # hberg_beech
  expect_equal(
    fisim:::xy_sample_regular(hberg_beech$boundary, n = 10, M = 1,
                              cell_size = NULL, random_rot = TRUE),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point = 1:11,
      x_s = c(119.19, 119.01, 102.11, 85.21, 68.31, 68.13, 51.23, 34.34, 34.15,
              17.26, 0.36),
      y_s = c(90.39, 5.54, 39.52, 73.49, 107.47, 22.62, 56.59, 90.57, 5.72,
              39.7, 73.67)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(hberg_beech$boundary, n = 10,  M = 1,
                              cell_size = NULL, random_rot = FALSE),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point = 1:9,
      x_s = c(26.04, 63.99, 101.94, 26.04, 63.99, 101.94, 26.04, 63.99, 101.94),
      y_s = c(13.11, 13.11, 13.11, 51.06, 51.06, 51.06, 89.01, 89.01, 89.01)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(hberg_beech$boundary, n = 10, M = 1,
                              cell_size = 1000, random_rot = FALSE),
    data.table::data.table(
      id_sample = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                    1L),
      id_point = 1:16,
      x_s = c(3.75, 35.38, 67, 98.62, 3.75, 35.38, 67, 98.62, 3.75, 35.38, 67,
              98.62, 3.75, 35.38, 67, 98.62),
      y_s = c(5.71, 5.71, 5.71, 5.71, 37.33, 37.33, 37.33, 37.33, 68.95, 68.95,
              68.95, 68.95, 100.58, 100.58, 100.58, 100.58)),
    tolerance = 5e-2)

  # kalimantan_peat
  expect_equal(
    fisim:::xy_sample_regular(kalimantan_peat$boundary, n = 10, M = 1,
                              cell_size = NULL, random_rot = TRUE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:9,
      x_s=c(72.47, 85.28, 98.1, 110.91, 35.75, 48.56, 61.38, 11.84, 24.66),
      y_s=c(3.5, 40.22, 76.94, 113.66, 16.32, 53.04, 89.76, 65.85, 102.57)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(kalimantan_peat$boundary, n = 10, M = 1,
                              cell_size = NULL, random_rot = FALSE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:9,
      x_s=c(38.37, 77.26, 116.15, 38.37, 77.26, 116.15, 38.37, 77.26, 116.15),
      y_s=c(36.57, 36.57, 36.57, 75.46, 75.46, 75.46, 114.35, 114.35, 114.35)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(kalimantan_peat$boundary, n = 10, M = 1,
                              cell_size = 2500, random_rot = FALSE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:6,
      x_s=c(11.06, 61.06, 111.06, 11.06, 61.06, 111.06),
      y_s=c(39.89, 39.89, 39.89, 89.89, 89.89, 89.89)),
    tolerance = 5e-2)

  # random_tree
  expect_equal(
    fisim:::xy_sample_regular(random_tree$boundary, n = 10, M = 1,
                              cell_size = NULL, random_rot = TRUE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:7,
      x_s=c(602.52, 665.38, 728.25, 387.93, 450.79, 173.34, 236.21),
      y_s=c(30.21, 244.79, 459.38, 93.07, 307.66, 155.94, 370.52)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(random_tree$boundary, n = 10, M = 1,
                              cell_size = NULL, random_rot = FALSE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:8,
      x_s=c(119.75, 343.35, 566.96, 790.57, 119.75, 343.35, 566.96, 790.57),
      y_s=c(220.86, 220.86, 220.86, 220.86, 444.46, 444.46, 444.46, 444.46)),
    tolerance = 5e-2)

  expect_equal(
    fisim:::xy_sample_regular(random_tree$boundary, n = 10, M = 1,
                              cell_size = 50000, random_rot = FALSE),
    data.table::data.table(
      id_sample=c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
      id_point=1:10,
      x_s=c(33.1, 256.71, 480.31, 703.92, 927.53, 33.1, 256.71, 480.31,
            703.92, 927.53),
      y_s=c(142.61, 142.61, 142.61, 142.61, 142.61, 366.22, 366.22, 366.22,
            366.22, 366.22)),
    tolerance = 5e-2)
})
