obj_to_mesh <- function(obj,
                          material = "gray") {
  obj_lines <- readLines(obj)
  
  vertex_lines <- obj_lines[grepl("^v ",obj_lines)]
  vertex_values <- as.numeric(unlist(strsplit(sub("v ","",vertex_lines)," ")))
  vertex_matrix <- t(matrix(vertex_values, nrow = length(vertex_lines), byrow = TRUE))
  vertex_matrix <- rbind(vertex_matrix, rep(1, ncol(vertex_matrix)))
  
  
  face_lines <- obj_lines[grepl("^f ",obj_lines)]
  face_values <- as.integer(sub("//.+","",unlist(strsplit(sub("f ","",face_lines)," "))))
  face_matrix <- t(matrix(face_values, nrow = length(face_lines), byrow = TRUE))
  
  mesh <- list(vb = vertex_matrix, it = face_matrix, primitivetype = "triangle")
  class(mesh) <- c("mesh3d", "shape3d")
  
  mesh
}
