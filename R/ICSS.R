##  ICSS -- Iterative Cumulative Sum of Squares algorithm by Inclan/Tiao (1994)
##
##  Copyright (C) 2020  Siegfried Köstlmeier <siegfried.koestlmeier@gmail.com>
##
##  This file is part of the ICSS R-package
##
##  ICSS is free software; you can redistribute it and/or
##  modify it under the terms of the GNU General Public License
##  as published by the Free Software Foundation; either version 2
##  of the License, or (at your option) any later version.
##
##  ICSS is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with ICSS If not, see <http://www.gnu.org/licenses/>.

#' Iterative Cumulative Sum of Squares (ICSS)
#' @description
#'    ICSS implements the Iterative Cumulative Sum of Squares (ICSS) algorithm by Inclan/Tiao (1994).
#'    The test detects structural breakpoints in the variance of time series data.
#'
#' @param data A numerical vector
#' @param demean An object of class "logical": If demean is TRUE, all data will get demeaned prior to the ICSS algorithm.
#'
#' @return ICSS returns a numerical vector containing the location of structural breakpoints or NA if none breakpoints are found.
#'
#' @export
ICSS <- function(data = data, demean = FALSE) {
  if (demean) data <- data - mean(data)

  S <- rstack::stack$new()
  S$push(c(1, length(data)))

  potential_change_points <- rep(0, length(data))
  cp_index <- 1

  while (!S$is_empty()) {
    current_range <- S$pop()
    print(current_range)
    change_points_sub <- ICSS_step_1_and_2(data[current_range[1]:current_range[2]])

    if (length(change_points_sub) == 2) {
      S$push(change_points_sub + current_range[1])
    }

    if (length(change_points_sub) > 0) {
      potential_change_points[cp_index:(cp_index + (length(change_points_sub) - 1))] <- change_points_sub + current_range[1]
    }

    cp_index <- cp_index + length(change_points_sub)
  }

  if (sum(potential_change_points, na.rm = T) == 0) warning("Unable to identify structural variance breaks in the series.")

  potential_change_points <- unique(c(0, sort(potential_change_points), length(data)))

  converged <- FALSE

  while (!converged) {
    new_cps_stack <- rstack::stack$new()
    for (i in 2:(length(potential_change_points) - 1)) {
      from <- potential_change_points[i - 1] + 1
      to <- potential_change_points[i + 1]

      if (is.na(to)) {
        return(NA)
      }

      Dk <- CenteredCusumValues(data[from:to])
      tmp <- check_critical_value(Dk)

      if (length(tmp$position) == 0) {
        return(NA)
      }
      exceeds <- ifelse(is.na(tmp$exceeds), FALSE, tmp$exceeds)
      position <- tmp$position

      if (exceeds) {
        new_cps_stack$push(from + position)
      }
    }

    stack_entries <- vector(mode = "numeric", length = new_cps_stack$size())
    for (j in 1:new_cps_stack$size()) {
      stack_entries[j] <- new_cps_stack$pop()
    }

    new_cps <- unique(c(0, sort(stack_entries), length(data)))
    converged <- is_converged(potential_change_points, new_cps)

    if (!converged) {
      potential_change_points <- new_cps
    }
  }

  change_points <- potential_change_points[2:(length(potential_change_points) - 1)]
  return(change_points)
}
