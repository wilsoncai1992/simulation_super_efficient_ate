indbegin <- which( colnames(data)=='hxb2.7.C.1mer')
indend <- which(colnames(data)=='hxb2.781.gap.1mer')	
AAcharactervars <- names(data)[indbegin:indend]
AAposns <- unlist (lapply (strsplit (AAcharactervars, split=".", fixed=T), function (x) x[2]))

VRC01contactsites <- 
c(97, 123, 124, 198, 276, 278, 279, 280, 281, 282, 365, 366, 367, 368, 371, 427, 428, 429, 430, 455, 456, 457, 458, 459, 460, 461, 463, 465, 466, 467, 469, 472, 473, 474, 476)
keepVRC01contactsites <- AAposns %in% VRC01contactsites
AAVRC01contactsitescharactervars <- AAcharactervars[keepVRC01contactsites]

make_majority_variant <- function(AAposn = 97, allAAposn = AAVRC01contactsitescharactervars){
	which_cols <- grep(paste0("hxb2.",AAposn), allAAposn)
	if(length(which_cols) == 0){
		return(NULL)
	}
	these_data <- data[, ..allAAposn][ , ..which_cols]
	max_var <- which.max(colSums(these_data))
	out <- these_data[,..max_var]
	colnames(out) <- paste0("hxb2.",AAposn,".majority")
	return(out)
}

all_majority <- sapply(VRC01contactsites, make_majority_variant, allAAposn = AAVRC01contactsitescharactervars)
null_boolean <- unlist(lapply(all_majority, is.null))
all_majority <- all_majority[!null_boolean]
all_cols <- Reduce(cbind, all_majority)

otherAAsites <- !(AAposns %in% VRC01contactsites)
othersites <- AAcharactervars[otherAAsites]
other_cols <- data[,..othersites]

data_rm_vrc01 <- data[,-..AAVRC01contactsitescharactervars]

Y <- as.numeric(data_rm_vrc01$ic50.censored)
W <- data_rm_vrc01[,-(1:3)]
A <- all_cols

make_data_set_for_analysis <- function(AAposn, allAAposn = AAVRC01contactsitescharactervars){
	A <- unlist(make_majority_variant(AAposn = AAposn, allAAposn = AAVRC01contactsitescharactervars), use.names = FALSE)
	Y <- as.numeric(data$ic50.censored)

	which_cols <- grep(paste0("hxb2.",AAposn), colnames(data))
	if(length(which_cols) == 0){
		return(NULL)
	}
	these_data <- data[ , -..which_cols][,-(1:3)]
	return(list(W = W, A = A, Y = Y))
}

all_data_sets_for_analysis <- sapply(VRC01contactsites, make_data_set_for_analysis)
names(all_data_sets_for_analysis) <- paste0("AAposn", VRC01contactsites)
null_data_sets <- unlist(lapply(all_data_sets_for_analysis, is.null))
all_data <- all_data_sets_for_analysis[!null_data_sets]

CD4bindingsites <- c(124, 125, 126, 127, 196, 198, 279, 280, 281, 282, 283, 365, 366, 367, 368, 369, 370, 374, 425, 426, 427, 428, 429, 430, 431, 432, 455, 456, 457, 458, 459, 460, 461, 469, 471, 472, 473, 474, 475, 476, 477)
keepCD4bindingsites <- AAposns %in% CD4bindingsites
AACD4bindingsitescharactervars <- AAcharactervars[keepCD4bindingsites]

