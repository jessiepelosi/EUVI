calculate.aban.un.stab.exp.pion <-
  function(grid.clim.pooled = grid.pooled, grid.clim.nat = grid.native, grid.clim.inv = grid.invaded, remove.marginal){
    
    
    # Remove marginal climatic conditions
    if(remove.marginal == 0){
      
      threshold.marginal.native <- 0
      threshold.marginal.invaded <- 0
      
    }else{
      
      threshold.marginal.native <- quantile(grid.clim.nat$Z[grid.clim.nat$Z > 0], probs = seq(0, 1, remove.marginal))[2]
      threshold.marginal.invaded <- quantile(grid.clim.inv$Z[grid.clim.inv$Z > 0], probs = seq(0, 1, remove.marginal))[2]
      
    }
    
    
    # Extract rasters per category
    abandonment_r<- (((grid.clim.nat$Z > threshold.marginal.native) - (grid.clim.inv$Z > threshold.marginal.invaded)) == 1) *
      (grid.clim.nat$z.uncor>0) *
      grid.clim.pooled$z.uncor
    pioneering_r <- (((grid.clim.inv$Z > threshold.marginal.invaded) - (grid.clim.nat$Z > threshold.marginal.native)) == 1) *
      (grid.clim.inv$z.uncor >0)*
      grid.clim.pooled$z.uncor 
    stability_r <- (((grid.clim.inv$Z > threshold.marginal.invaded) * (grid.clim.nat$Z > threshold.marginal.native)) == 1) * 
      ((grid.clim.inv$z.uncor * grid.clim.nat$z.uncor)>0)*
      grid.clim.pooled$z.uncor
    expansion_r <- (((grid.clim.inv$Z > threshold.marginal.invaded) * (grid.clim.nat$Z > threshold.marginal.native)) == 1) *
      ((grid.clim.inv$z.uncor >0) * (grid.clim.nat$z.uncor==0))*
      grid.clim.pooled$z.uncor 
    unfilling_r<-(((grid.clim.inv$Z > threshold.marginal.invaded) * (grid.clim.nat$Z > threshold.marginal.native)) == 1) *
      ((grid.clim.nat$z.uncor >0) * (grid.clim.inv$z.uncor==0))*
      grid.clim.pooled$z.uncor 
    
    # Calculate proportions
    abandonment <- terra::global(abandonment_r, "sum") / terra::global(grid.clim.pooled$z.uncor, "sum")
    pioneering <- terra::global(pioneering_r, "sum") / terra::global(grid.clim.pooled$z.uncor, "sum")
    stability<- terra::global(stability_r, "sum") / terra::global(grid.clim.pooled$z.uncor, "sum")
    expansion<- terra::global(expansion_r, "sum") / terra::global(grid.clim.pooled$z.uncor, "sum")
    unfilling<- terra::global(unfilling_r, "sum") / terra::global(grid.clim.pooled$z.uncor, "sum")
    
    result <- c(abandonment, unfilling, stability, expansion, pioneering)
    names(result) <- c("abandonment", "unfilling", "stability", "expansion", "pioneering")
    
    return(result)
    
    
  }