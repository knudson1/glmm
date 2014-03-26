getFamily <-
function(family.mcml)
{
	if(is.character(family.mcml))
		family.mcml<-get(family.mcml,mode="function",envir=parent.frame())
	if(is.function(family.mcml))
		family.mcml<-family.mcml()
	#right way to check the class if there might be more than one class (maybe bc of hierarchies)
	if(!inherits(family.mcml,"mcml.family")) 
		stop(" 'family' not recognized") 
	
	family.mcml	
}
