################
# mdf_par_name #
################

mdf_par_name <- function(nm, pb_chr = c('-', ' ', '/', '*', '+', '_')){
  
  for(c in 1:length(pb_chr)){nm <- gsub(pattern = pb_chr[c],
                            replacement = "", x = nm)}
  return(nm)
}