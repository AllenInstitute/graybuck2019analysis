library(rgl)
library(cocoframer)
library(purrr)
options(stringsAsFactors = F)
source("read_aibs_obj.R")

obj_to_mesh <- function(obj,
                        material = "gray",
                        invert_y = TRUE,
                        yrange = c(0,1)) {
  
  obj_lines <- readLines(obj)
  
  vertex_lines <- obj_lines[grepl("^v ",obj_lines)]
  vertex_values <- as.numeric(unlist(strsplit(sub("v ","",vertex_lines)," ")))
  vertex_matrix <- t(matrix(vertex_values, nrow = length(vertex_lines), byrow = TRUE))
  
  if(invert_y) {
    midpoint <- (yrange[1] + yrange[2]) / 2
    vertex_matrix[2,] <- vertex_matrix[2,] + 2 * (midpoint - vertex_matrix[2,])
  }
  
  vertex_matrix <- rbind(vertex_matrix, rep(1, ncol(vertex_matrix)))
  
  
  face_lines <- obj_lines[grepl("^f ",obj_lines)]
  face_values <- as.integer(sub("//.+","",unlist(strsplit(sub("f ","",face_lines)," "))))
  face_matrix <- t(matrix(face_values, nrow = length(face_lines), byrow = TRUE))
  
  mesh <- list(vb = vertex_matrix, it = face_matrix, primitivetype = "triangle", material = material)
  class(mesh) <- c("mesh3d", "shape3d")
  
  mesh
}

obj_dir <- "http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_meshes/"

ont <- get_mba_ontology()
ont <- flatten_mba_ontology(ont)

structures <- c("root","grey","CH","CB","Isocortex","VISp","VISp5","DORpm","sAMY","MOs","APr","P-mot","CBN","BLA","SSp","RHP","HATA","cc","MOp","CLA")

structure_ids <- ont$id[match(structures, ont$acronym)]
structure_urls <- paste0(obj_dir, structure_ids, ".obj")
structure_files <- paste0(structure_ids, ".obj")

walk(structure_urls,
     function(x) {
       if(!file.exists(sub(".+/","",x))) {
         download.file(x,
                       sub(".+/","",x))
       }
       
     })

structure_meshes <- map(structure_files, 
                        obj_to_mesh,
                        invert_y = TRUE,
                        yrange = c(140.027, 7558.69))

names(structure_meshes) <- structures

plot_structure <- function(mesh_list,
                           main_structure,
                           main_color = "#74CAFF",
                           main_alpha = 1,
                           background_structure = NULL,
                           background_color = "#808080",
                           background_alpha = 0.2) {
  if(is.null(background_structure)) {
    meshes <- mesh_list[main_structure]
  } else {
    meshes <- mesh_list[c(main_structure,
                          background_structure)]
    
    meshes[[background_structure]]$material <- list(color = background_color,
                                                    alpha = background_alpha)
  }
  
  meshes[[main_structure]]$material <- list(color = main_color,
                                                  alpha = main_alpha)
  
  view3d(theta = -45, phi = 35, zoom = 0.7)
  
  shapelist3d(meshes,
              box = FALSE,
              axes = FALSE,
              xlab = "", ylab = "", zlab = "")
}

open3d()
walk(structures,
     function(x) {
       if(!file.exists(paste0(x,".png"))) {
         clear3d()
         
         if(x == "root") {
           plot_structure(structure_meshes,
                          main_structure = x)
         } else {
           plot_structure(structure_meshes,
                          main_structure = x,
                          background_structure = "root")
         }
         rgl.snapshot(paste0(x,".png"))
         
       }
       
     })

plot_structure(structure_meshes,
               main_structure = "grey")

rgl.snapshot("grey.png")

clear3d()
plot_structure(structure_meshes,
               main_structure = "VISp",
               background_structure = "grey")

rgl.snapshot("VISp.png")

clear3d()
plot_structure(structure_meshes,
               main_structure = "VISp",
               background_structure = "grey")

rgl.snapshot("VISp.png")



shapelist3d(structure_meshes[c("grey","VISp")], 
       alpha = 0.2,
       box = FALSE,
       axes = FALSE,
       xlab = "", ylab = "", zlab = "")

sp <- spin3d(axis = c(0.5, 0.5, 0),
             rpm = 2)

play3d(spin3d(rpm = 2), duration = 2)

pars <- par3d()


prj <- par3d("userMatrix") 
prjs <- rowSums(prj) 
theta <- acos(prjs[3]/sqrt(sum(prjs[1:3]^2))) * 180 / pi 
phi <- atan(prjs[2]/prjs[1]) * 180/pi
