
#' plot_barcode
#'
#' Plotting barcode using the TDAstats function.
#'
#' @param barcode_individual A barcode matrix derived from the compute_homology or compute_homology_multi function.
#'
#' @return a plot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_dimensions(barcode_individual)
#' }
plot_barcodes <- function(barcode_individual){
  return(TDAstats::plot_barcode(barcode_individual))
}


#' plot_dimensions
#'
#' Function that plots the barcode and the PBFs for a specific barcode.
#'
#' @param barcode_df Barcode to be plotted
#' @param resolution Resolution, interval width, for the construction of the Persistent Betty Function.
#'
#' @return Returns a plot of the barcodes and the Persistent Betti Functions.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_dimensions(barcode_df,0.1)
#' }
plot_dimensions <- function(barcode_df,resolution = 0.1){
  PBF_barcode_zero <- compute_pbf_for_curve(seq(from = 0, to = max(barcode_df[,2:3])+5,by = resolution),barcode_df,k = 2,dim_betti = 0)
  df_bc_zero <- base::data.frame(x = seq(from = 0, to = max(barcode_df[,2:3])+5,by = resolution),y =PBF_barcode_zero)
  p_zero <- ggplot2::ggplot(df_bc_zero, ggplot2::aes(x=x, y=y)) + ggplot2::geom_line(color="red") + ggplot2::theme_bw() +  ggplot2::theme(axis.title.x=ggplot2::element_blank()) +  ggplot2::ylab(expression(beta[0]))
  PBF_barcode_one <- compute_pbf_for_curve(seq(from = 0, to = max(barcode_df[,2:3])+5,by = resolution),barcode_df,k = 2,dim_betti = 1)
  df_bc_one <- base::data.frame(x = seq(from = 0, to = max(barcode_df[,2:3])+5,by =resolution),y =PBF_barcode_one)
  p_one <- ggplot2::ggplot(df_bc_one, ggplot2::aes(x=x, y=y)) + ggplot2::geom_line(color="green") + ggplot2::theme_bw() + ggplot2::theme(axis.title.x=ggplot2::element_blank()) +  ggplot2::ylab(expression(beta[1]))
  PBF_barcode_two <- compute_pbf_for_curve(seq(from = 0, to = max(barcode_df[,2:3])+5,by = resolution),barcode_df,k = 2,dim_betti = 2)
  df_bc_two <- base::data.frame(x = seq(from = 0, to = max(barcode_df[,2:3])+5,by = resolution),y =PBF_barcode_two)
  p_two <- ggplot2::ggplot(df_bc_two, ggplot2::aes(x=x, y=y)) + ggplot2::geom_line(color="blue") + ggplot2::theme_bw() + ggplot2::theme(axis.title.x=ggplot2::element_blank()) +  ggplot2::ylab(expression(beta[2]))
  plot_arranged <- cowplot::plot_grid(p_two, p_one, p_zero,ncol = 1, nrow = 3)
  return(plot_arranged)
}
