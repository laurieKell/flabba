
validity<-function(object) {
  return(TRUE)
  ## Catch must be continous
  yrs<-dimnames(catch(object))$year
  
  if (!all(yrs == ac(dims(catch(object))$minyear:dims(catch(object))$maxyear)))
    return("years in catch not continous")
  
  # range
  dims <-dims(object)
  range<-as.list(object@range)
  
  return(TRUE)}

#' @title jabba class
#'
#' @description A class that implement a biomass dynamic stock assessment model using JABBA. 
#' 
#' @details 
#' ...
#'  
#' @slot name     {A \code{character} name of  stock}
#' @slot desc     {A \code{character} providing a fuller description of the object}       
#' @slot version  {A \code{character} version}
#' @slot range    {A \code{numeric} vector containing the quant and year ranges}
#' @slot catch    {An \code{FLQuant}  total catch by year}        
#' @slot index    {An \code{FLQuants} indicies of abundance}        
#' @slot index.se {An \code{FLQuants} SEs for indicies of abundance}        
#' @slot stock    {An \code{FLQuant} estimated stock by year}       
#' @slot priors   {An \code{array} }       
#' @slot diags    {A \code{data.frame} with residuals and covariates from fit of CPUE to stock }     
#' 
#' All slots in the class have accessor and replacement methods that provide validation and protection of their data.
#' 
#' @export
#' @import FLCore 
#' @import methods
#'
#' @aliases 
#' 
#' @rdname jabbaClass
#'  
#' @examples
#' \dontrun{jabba()}
 setClass('biodyn', representation('FLComp',
    version       ='character',
    catch         ='FLQuant',  
    indices       ='FLQuants',
    indices.se    ='FLQuants',
    stock         ='FLQuant',
    diags         ='data.frame',
    priors        ='array'
    ),
  prototype(
    version     ="v1.6beta",
    range       =unlist(list(minyear=as.numeric(NA), maxyear=as.numeric(NA))),
    catch       =FLQuant(),
    stock       =FLQuant(),
    indices     =FLQuants(),
    indices.se  =FLQuants(),
    diags       =data.frame()), 
	validity=validity)
