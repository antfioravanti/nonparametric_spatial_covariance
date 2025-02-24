# Extract reference coordinates
  ref_x = point_df$x[1]
  ref_y = point_df$y[1]
  distances = sqrt((grid$x - ref_x)^2 + (grid$y - ref_y)^2)
  same_point_indices = which(grid$x == ref_x & grid$y == ref_y)
  grid$dist = distances
  # Identify the Closest Point
  closest_index = order(distances)[2]
  closest_point = grid[closest_index, ]
  
  # Identify the Medium-Distance Point
  median_distance = median(distances)
  median_diff = abs(distances - median_distance)
  medium_index = which.min(median_diff)
  medium_point = grid[medium_index, ]
  
  # Identify the Farthest Point
  farthest_index = which.max(distances)
  farthest_point = grid[farthest_index, ]
  